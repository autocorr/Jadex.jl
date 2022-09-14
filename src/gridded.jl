module GridRunner

export runseq, rungrid, get_interp, interp_errors

using Base.Threads
using DataFrames
using Interpolations

using ..EscapeProbability
using ..ReadData: Specie
using ..Background
using ..RunDefinition
using ..Solver


function draw_from_axes(axes)
    N = length(axes)
    draws = Array{Float64}(undef, N)
    for i in eachindex(axes)
        vmin, vmax = extrema(axes[i])
        draws[i] = vmin + (vmax - vmin) * rand()
    end
    draws
end


"""
    runseq(sol, rdfs; ...)

Execute a series model definitions. All `RunDef` run definitions should be
for the same `Specie` type (i.e., LAMDA file).
"""
function runseq(sol::Solution, rdfs::AbstractVector{<:RunDef};
        freq_min=-Inf, freq_max=Inf)
    @assert length(rdfs) >= nthreads()
    mol_name = rdfs[1].mol.name
    sols = [deepcopy(sol) for _ in 1:nthreads()]
    dfs = DataFrame[]
    lk = ReentrantLock()
    @threads for i in eachindex(rdfs)
        rdf = rdfs[i]
        sol_t = sols[threadid()]
        # Ensure that all `RunDef`s are for the same `Specie`, and thus
        # have `Solution` instances with consistent array dimensions.
        @assert rdf.mol.name == mol_name
        solve!(sol_t, rdf)
        df = get_results(sol_t, rdf; freq_min=freq_min, freq_max=freq_max)
        # Ensure that there is not a data race in appending the results.
        lock(lk) do
            push!(dfs, df)
        end
        reset!(sol_t)
    end
    vcat(dfs...; cols=:union)
end

function runseq(rdfs::AbstractVector{<:RunDef}; kwargs...)
    sol = Solution(rdfs[1])
    runseq(sol, rdfs; kwargs...)
end


"""
    rungrid(mol, density, tkin, cdmol, deltav, transitions; ...)

Generate a grid of τ and Tₑₓ values over the set of physical conditions and
selected transitions. Most parameters for a `RunDef` run definition can be
set as keyword arguments (e.g., β, bg). Arguments can be passed as scalars
or collections. Runs will be executed in parallel if Julia was started
with multiple threads.

Resulting cubes are indexed with shape `(n, Tₖ, N, δv, j)` for volume density,
kinetic temperature, column density, linewidth, and transition number. The
transition index number corresponds to the TRANS value in the LAMDA data file
or the index number in the output table.
"""
function rungrid(mol::Specie,
        density::AbstractVector,
        tkin::AbstractVector,
        cdmol::AbstractVector,
        deltav::AbstractVector,
        transitions::AbstractVector{I};
        collision_partner::String="h2",
        escprob::Function=βlvg,
        bg::Union{Nothing, BackgroundField}=nothing,
        reduced::Bool=false,
        minpop=1e-20,
        tolerance=1e-6,
        iterpars=IterationParams(),
    ) where {I <: Integer}
    # Validate radiative transition indices
    valid_transitions = collect(1:mol.transitions.n)
    @assert all(transitions .∈ Ref(valid_transitions))
    # Validate collision partner names
    # FIXME Should add a special case for (h2)->(p-h2,o-h2)
    #       for certain species using `thermal_h2_density`.
    valid_partners = [c.name for c in mol.colliders]
    @assert all(collision_partner .∈ Ref(valid_partners))
    # Pre-allocate background if given as nothing
    bg = bg === nothing ? galactic_isrf(mol) : bg
    # Create output data cubes
    param_shape = length.((density, tkin, cdmol, deltav))
    cube_shape = (param_shape..., length(transitions))
    tau_cube = zeros(cube_shape)
    tex_cube = zeros(cube_shape)
    # Compute values over cube
    @assert prod(param_shape) >= nthreads()
    sols = [Solution(mol, iterpars) for _ in 1:nthreads()]
    @threads for indices in CartesianIndices(param_shape)
        i_n, i_T, i_N, i_δ = Tuple(indices)
        n = density[i_n]
        T = tkin[i_T]
        N = cdmol[i_N]
        δ = deltav[i_δ]
        n_dict = Dict(collision_partner => n)
        rdf = RunDef(mol, density=n_dict, tkin=T, cdmol=N, deltav=δ,
                     escprob=escprob, bg=bg, reduced=reduced, minpop=minpop,
                     tolerance=tolerance)
        sol = sols[threadid()]
        solve!(sol, rdf)
        for (i, j) in enumerate(transitions)
            tau_cube[i_n, i_T, i_N, i_δ, i] = sol.τl[j]
            tex_cube[i_n, i_T, i_N, i_δ, i] = sol.tex[j]
        end
        reset!(sol)
    end
    tau_cube, tex_cube
end

function rungrid(mol::Specie, args...; kwargs...)
    args = [arg isa Number ? [arg] : collect(arg) for arg in args]
    rungrid(mol, args...; kwargs...)
end
rungrid(mol::Specie, args; kwargs...) = rungrid(mol, args...; kwargs...)


"""
    get_interp(A::Matrix, axes; B=<interpolation-type>)

Generate a 5D interpolation object for a given τ or Tₑₓ cube where the last
dimension indexes the transition number. By default a second order (quadratic)
interpolation is used with flat extrapolation for out-of-bounds access. The
axis containing the transition number is not interpolated -- an exactly
matching value must be given.

Please ensure that the axes passed have a linear step size (nb. use log(n) or
log(N) axes if using logarithmically spaced values). A gridded interpolation
object can be created using the `Gridded` interpolation type.
"""
function get_interp(A::AbstractArray{<:Real}, axes;
        B=BSpline(Quadratic(Line(OnGrid()))))
    # `scale` requires that the axes be given as ranges, not vectors
    axis_ranges = [
            range(start=ax[1], stop=ax[end], length=length(ax))
            for ax in axes[1:end-1]
    ]
    push!(axis_ranges, UnitRange(1, length(axes[end])))
    # Only look-up for transition axis
    interp_types = (repeat([B], length(axes)-1)..., NoInterp())
    itp = interpolate(A, interp_types)
    scale(itp, axis_ranges...)
end


"""
Compare the interpolation to a uniform random sample and return an array
of the deviations.

Note that this implementation currently only works for 5 parameter grids
of (n, Tₖ, N, Δv, j) and where `n` is the volume density of the specified
collision partner.
"""
function interp_errors(itp, axes, mol, trans;
        property=:τl, inlog=nothing, nsamples=10_000, collision_partner="h2",
        kwargs...)
    # FIXME `trans` refers to different values in `itp` and index in the property.
    param_axes = axes[1:end-1]
    transitions = axes[end]
    trans_index = findfirst(==(trans), transitions)
    nparams = length(param_axes)
    log_by_val = inlog === nothing ? repeat([false], nparams) : inlog
    deviations = Float64[]
    lk = ReentrantLock()
    @threads for _ in 1:nsamples
        params = draw_from_axes(param_axes)
        interp_v = itp(params..., trans_index)
        # Compute precise value of grid point
        lin_params = [l ? exp10(p) : p for (l, p) in zip(log_by_val, params)]
        density = Dict(collision_partner => lin_params[1])
        rdf = RunDef(mol; density=density, tkin=lin_params[2], cdmol=lin_params[3],
                     deltav=lin_params[4], kwargs...)
        _, sol = solve(rdf)
        true_v = getproperty(sol, property)[trans]
        rel_diff = abs(true_v - interp_v) / true_v
        lock(lk) do
            push!(deviations, rel_diff)
        end
    end
    deviations
end


end  # module
