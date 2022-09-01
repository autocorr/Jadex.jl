module GridRunner

export rungrid

using Base.Threads
using Interpolations

using ..EscapeProbability
using ..ReadData: Specie
using ..Background
using ..RunDefinition
using ..Solver


function runseq(sol::Solution, rdfs::AbstractVector{RunDef};
        min_freq=-Inf, max_freq=Inf)
    error("Not implemented.")
    # FIXME
    # - add multi-threading
    # NOTE
    # - solutions are `Specie` dependent, so have to be careful when multiple
    #   molecules are run from the same list. This should be caught by
    #   `validate`, so hopefully will avoid a bounds error, but still needs to
    #   be handled carefully. Should not expose this interface as the primary
    #   one for this reason.
    dfs = DataFrame[]
    for rdf in rdfs
        # FIXME Check if `rdf` has a different molecule type, if so, then clone
        # a new solution instance so the same parameters are kept.
        solve!(sol, rdf)
        push!(dfs, get_results(sol, rdf; min_freq=min_freq, max_freq=max_freq))
        reset!(sol)
    end
    vcat(dfs...; cols=:union)
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
    @assert prod(param_shape) > nthreads()
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
        for j in transitions
            tau_cube[i_n, i_T, i_N, i_δ, j] = sol.τl[j]
            tex_cube[i_n, i_T, i_N, i_δ, j] = sol.tex[j]
        end
        reset!(sol)
    end
    tau_cube, tex_cube
end

function rungrid(mol::Specie, args...; kwargs...)
    args = [arg isa Number ? [arg] : collect(arg) for arg in args]
    rungrid(mol, args...; kwargs...)
end


"""
    create_interp(A::Matrix, axes; B=<interpolation-type>)

Generate a 5D interpolation object for a given τ or Tₑₓ cube where the last
dimension indexes the transition number. By default a second order (quadratic)
interpolation is used with flat extrapolation for out-of-bounds access. The
axis containing the transition number is not interpolated -- an exactly
matching value must be given.

Please ensure that the axes passed have a linear step size (nb. use log(n) or
log(N) axes if using logarithmically spaced values). A gridded interpolation
object can be created using the `Gridded` interpolation type.
"""
function create_interp(A::Matrix, axes; B=BSpline(Quadratic(Flat(OnGrid()))))
    # `scale` requires that the axes be given as ranges, not vectors
    axis_ranges = [range(ax[1], ax[end], step=ax[2]-ax[1]) for ax in axes]
    # Only look-up for transition axis
    itp = interpolate(axes, A, (B, B, B, B, NoInterp()))
    scale(itp, axis_ranges)
end


end  # module
