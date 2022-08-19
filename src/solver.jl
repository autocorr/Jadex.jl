module Solver

export Solution
export reset!, solve!, get_results

using Base.Threads
using Logging
using LinearAlgebra
using DataFrames

using ..Constants: FK, THC, CLIGHT, KBOLTZ
using ..RunDefinition: RunDef


const GAUSS_FWHM = √(8 * log(2))        #  ~2.355
const GAUSS_AREA = √(2π) / GAUSS_FWHM   #  ~1.065
const FGAUSS = GAUSS_AREA * 8π          # ~26.753
# Below is the precise Float64 value used by RADEX, which is accurate to about
# six significant figures given the approximation `√π/2√log2 ~ 1.0645`.
#const FGAUSS = √π / (2 * √log(2)) * 8π
#const FGAUSS = 26.753802360251857


struct IterationParams
    min::Int
    max::Int
    start::Int
    pending::Int

    function IterationParams(min, max, start, pending)
        @assert 0 <= min < max
        @assert start >= 4
        @assert pending >= 4
        new(min, max, start, pending)
    end
end
IterationParams() = IterationParams(5, 10_000, 4, 4)


mutable struct Solution{F <: AbstractFloat}
    iter::IterationParams
    niter::Int
    nthick::Int
    τsum::F
    rhs::Vector{F}
    yrate::Matrix{F}
    tex::Vector{F}
    xpop::Vector{F}
    xpop_::Matrix{F}
    τl::Vector{F}
end

function Solution(rdf::RunDef{F,B}, iter::IterationParams) where {F, B}
    niter = 0
    nthick = 0
    τsum = zero(F)
    mol = rdf.mol
    nlev = mol.levels.n
    ntran = mol.transitions.n
    # Initialize Y rate matrix and right-hand-side vector
    rhs = zeros(F, nlev)
    rhs[end] = one(eltype(rhs))  # for conservation eq.
    yrate = zeros(F, nlev, nlev)
    clear_rates!(yrate)
    # Initialize excitation temperature and population level arrays
    τl = zeros(F, ntran)
    tex = zeros(F, ntran)
    xpop = zeros(F, nlev)
    xpop_ = zeros(F, (nlev, 3))
    # Default constructor for Solution type
    Solution(iter, niter, nthick, τsum, rhs, yrate, tex, xpop, xpop_, τl)
end
Solution(rdf::RunDef) = Solution(rdf, IterationParams())


function clear_rates!(yrate)
    yrate .= zero(eltype(yrate))
end
clear_rates!(sol::Solution) = clear_rates!(sol.yrate)


function reset!(sol::Solution; keep_result=false)
    sol.niter = zero(sol.niter)
    sol.nthick = zero(sol.nthick)
    sol.τsum = zero(sol.τsum)
    clear_rates!(sol.yrate)
    sol.xpop_ .= zero(eltype(sol.xpop_))
    if !keep_result
        sol.τl .= zero(eltype(sol.τl))
        sol.tex .= zero(eltype(sol.tex))
        sol.xpop .= zero(eltype(sol.xpop))
    end
    sol
end


function validate(sol::Solution, rdf::RunDef)
    nlev = rdf.mol.levels.n
    ntran = rdf.mol.transitions.n
    @assert length(sol.τl) == ntran
    @assert length(sol.tex) == ntran
    @assert length(sol.rhs) == nlev
    @assert length(sol.xpop) == nlev
    @assert size(sol.xpop_) == (nlev,    3)
    @assert size(sol.yrate) == (nlev, nlev)
    return nothing
end


"""
On the first iteration, use the background intensity to estimate the
level populations from optically thin statistical equillibrium. This
assumes that the emission is optically thin by setting the escape
probability β=1.
"""
function init_radiative!(sol::Solution, rdf::RunDef)
    l = rdf.mol.levels
    t = rdf.mol.transitions
    yrate = sol.yrate
    @inbounds for (i, m, n) in zip(1:t.n, t.iupp, t.ilow)
        gmn = l.gstat[m] / l.gstat[n]
        etr = FK * t.xnu[i] / rdf.bg.trj[i]
        exr = etr >= 160.0 ? zero(etr) : (inv ∘ expm1)(etr)
        A = t.aeinst[i]
        yrate[m,m] += A * (1 + exr)
        yrate[n,n] += A * gmn * exr
        yrate[m,n] -= A * gmn * exr
        yrate[n,m] -= A * (1 + exr)
    end
    sol
end


"""
    step_radiative!(sol::Solution, rdf::RunDef)

Compute contribution of radiative processes to the rate matrix using the
escape probability. Modifies the `Solution` in-place.
"""
function step_radiative!(sol::Solution, rdf::RunDef)
    l = rdf.mol.levels
    t = rdf.mol.transitions
    yrate = sol.yrate
    cddv = rdf.cdmol / rdf.deltav
    xpop = sol.xpop
    τl   = sol.τl
    # Count optically thick lines
    nthick = 0
    @inbounds for (i, xnu, m, n, A, tb) in zip(
                1:t.n, t.xnu, t.iupp, t.ilow, t.aeinst, rdf.bg.totalb)
        xt = xnu^3
        gmn = l.gstat[m] / l.gstat[n]
        # Calculate line optical depth
        pop_factor = xpop[n] * gmn - xpop[m]
        τl[i] = cddv * pop_factor / (FGAUSS * xt / A)
        if τl[i] > 1e-2; nthick += 1 end
        # Use escape probability approximation for internal intensity.
        # Note that excluding the contribution to the local radiation
        # field is effectively lambda acceleration for a single-zone
        # model.
        β   = rdf.escprob(τl[i])
        exr = tb * β / (THC * xt)
        # Radiative contribution to the rate matrix
        yrate[m,m] += A * (β + exr)
        yrate[n,n] += A * gmn * exr
        yrate[m,n] -= A * gmn * exr
        yrate[n,m] -= A * (β + exr)
    end
    sol.nthick = nthick
    sol
end


"""
    step_collision!(sol::Solution, rdf::RunDef)

Compute contribution of collisional processes to the rate matrix. Modifies
the `Solution` in-place.
"""
function step_collision!(sol::Solution, rdf::RunDef)
    nlev = rdf.mol.levels.n
    yrate = sol.yrate
    crate = rdf.crate
    ctot  = rdf.ctot
    # About ~10% faster to loop (j,i) than (i,j) -- `setindex!` more
    # performance sensitive than `getindex`.
    @inbounds for j = 1:nlev
        yrate[j,j] += ctot[j]
        for i = 1:nlev
            if i != j
                yrate[i,j] -= crate[j,i]
            end
        end
    end
    sol
end


@inline function store_population!(xpop_, xpop)
    xpop_[:, 3] .= @view xpop_[:, 2]
    xpop_[:, 2] .= @view xpop_[:, 1]
    xpop_[:, 1] .= xpop
    xpop_
end
@inline function store_population!(sol::Solution)
    store_population!(sol.xpop_, sol.xpop)
    sol
end


@doc raw"""
    solve_rate_eq!(xpop, yrate, rhs)

Solve the system of linear equations for the excitation rates to and from each
energy level. The result is computed using LU factorization and back-substition
on the `yrate` matrix without pivoting.

The system of equations ``Y x = r`` is inverted to ``x = Y^{-1} r`` and solved,
where ``Y`` (`yrate`) is the rate matrix encoding the total rate of
excitation/dexcitation from one level to another, ``x`` (`xpop`) is the level
populations, and ``r`` is the right-hand side (`rhs`) expressing the time
rate-of-change of the level populations. Statistical equillibrium assumes that
the derivative `rhs` term is all zero except for the conservation equation.

## Notes
These operations correspond to the LU factorization subroutines `ludcmp`
(lower-upper decomposition) and `lubksb` (lower-upper back substition) Fortran
functions in the F77 RADEX source code `matrix.f` and `slatec.f`. In
practice these subroutines simply prepare (and copy) the arguments for the
subroutine `SGEIR`, which does both the decomposition and solution through
back-substitution.  In `lubksb`, the equation for the last energy level is
over-written with the conservation equation (i.e., ones) to make the rate
matrix square and non-singular (otherwise the last "column" would be all
zeros). This is also done here by over-writing the rates in-place with
ones.

The LU factorization `lu!` replaces `yrate` in-place and returns a view.
Some work products are still generated with LAPACK/BLAS, but benchmarking
shows that these impose negligible overhead.

A more rigorous factorization would be to include the conservation equation as
an additional row and solve the rectangular system of equations ``(m+1,m)``.
This is also what `DESPOTIC` does using the non-negative least squares (nnls)
solver. Tests with QR factorization and Gaussian elimination however were less
numerically stable than LU factorization as well as a factor of a few slower.
"""
function solve_rate_eq!(xpop, yrate, rhs)
    # Add conservation equation
    yrate[end,:] .= one(eltype(yrate))
    # These operations are performed in-place and are equivalent to:
    #   xpop .= yrate \ rhs
    Q = try
        lu!(yrate, NoPivot())
    catch e
        if isa(e, ZeroPivotException)
            lu!(yrate)
        else
            rethrow()
        end
    end
    ldiv!(xpop, Q, rhs)
end
function solve_rate_eq!(sol::Solution)
    solve_rate_eq!(sol.xpop, sol.yrate, sol.rhs)
    sol
end


@doc raw"""
    solve_rate_eq_reduced!(xpop, yrate, nr)

Solve the rate equation for a restricted/reduced number of levels. Levels
higher than `nr` are assumed to be optically thin and only coupled by radiative
processes. A cascade to levels less than or equal to `nr` is computed for
higher-lying transitions.

This formalism may not be strictly accurate for complex molecules (e.g.,
``\mathrm{CH_3OH}``) where many transitions between levels may be radiatively
forbidden.
"""
function solve_rate_eq_reduced!(xpop, yrate, nr)
    Y = yrate
    nlev = size(Y, 2)
    dtype = eltype(Y)
    @assert nr > 1
    @assert nr + 1 < nlev
    # Create reduced matrix and add normalized contributions from outside
    # levels.
    Yr = Y[1:nr+1,1:nr+1]  # slicing copies
    @inbounds for j = 1:nr, i = 1:nr, k = nr+1:nlev
        Yr[i,j] += abs(Y[k,j] * Y[i,k] / Y[k,k])
    end
    # Allocate or create views for reduced RHS and level populations.
    # Note that the (n+1) level population gets overwritten by the
    # conservation equation, but gets put back in when computing the cascade.
    rhs = zeros(dtype, nr+1)
    rhs[end] = one(dtype)
    x = @view xpop[1:nr+1]
    solve_rate_eq!(x, Yr, rhs)
    # Compute cascade for higher excitation lines. This is not accessed
    # efficiently, but it appears fixed by the problem (total population
    # from all higher states to a given lower state).
    @inbounds for k = nr+1:nlev
        xtot = zero(dtype)
        for j = 1:nr
            xtot += xpop[j] * Y[k,j]
        end
        xpop[k] = abs(xtot / Y[k,k])
    end
    xpop
end
function solve_rate_eq_reduced!(sol::Solution, rdf::RunDef)
    solve_rate_eq_reduced!(sol.xpop, sol.yrate, rdf.nreduced)
    sol
end


@inline function floor_population!(xpop, minpop::Real)
    @. xpop = max(minpop, xpop)
end
@inline function floor_population!(sol::Solution, rdf::RunDef)
    floor_population!(sol.xpop, rdf.minpop)
    sol
end


@inline function islowpop(xpop, n::Integer, m::Integer, minpop::Real)
    xpop[n] <= minpop || xpop[m] <= minpop
end


function init_tau_tex!(sol::Solution, rdf::RunDef)
    t = rdf.mol.transitions
    tex   = sol.tex
    xpop  = sol.xpop
    gstat = rdf.mol.levels.gstat
    @inbounds for (i, xnu, m, n) in zip(1:t.n, t.xnu, t.iupp, t.ilow)
        if islowpop(xpop, n, m, rdf.minpop)
            tex[i] = rdf.bg.totalb[i]
        else
            pop_factor = xpop[n] * gstat[m] / (xpop[m] * gstat[n])
            tex[i] = FK * xnu / log(pop_factor)
        end
    end
    sol
end


"""
    step_tau_tex!(sol::Solution, rdf::RunDef)

Compute the optical depth and excitation temperature for each line.
"""
function step_tau_tex!(sol::Solution, rdf::RunDef)
    t    = rdf.mol.transitions
    τl   = sol.τl
    tex  = sol.tex
    xpop = sol.xpop
    gstat = rdf.mol.levels.gstat
    cddv = rdf.cdmol / rdf.deltav
    τsum = 0.0
    @inbounds for (i, xnu, m, n) in zip(1:t.n, t.xnu, t.iupp, t.ilow)
        if islowpop(xpop, n, m, rdf.minpop)
            itex = tex[i]
        else
            pop_factor = xpop[n] * gstat[m] / (xpop[m] * gstat[n])
            itex = FK * xnu / log(pop_factor)
        end
        # Only optically thick lines count towards convergence
        if τl[i] > 0.01
            τsum += abs((itex - tex[i]) / itex)
        end
        # Update excitation temperature and optical depth
        tex[i] = 0.5 * (itex + tex[i])
        pop_factor = xpop[n] * gstat[m] / gstat[n] - xpop[m]
        τl[i] = cddv * pop_factor / (FGAUSS * xnu^3 / t.aeinst[i])
    end
    sol.τsum = τsum
    sol
end


"""
    ng_accelerate!(sol::Solution)

Apply Ng-acceleration to level populations from previous three iterations.  See
S4.4.7 of "Radiative Transfer in Astrophysics: Theory, Numerical Methods, and
Applications" lecture note series by C.P. Dullemond which is itself based on
Olson, Auer, and Buchler (1985).
"""
function ng_accelerate!(sol::Solution)
    x1 = sol.xpop
    x2 = @view sol.xpop_[:, 1]
    x3 = @view sol.xpop_[:, 2]
    x4 = @view sol.xpop_[:, 3]
    q1 = @. x1 - 2x2 + x3
    q2 = @. x1 -  x2 - x3 + x4
    q3 = @. x1 -  x2
    A1 = q1 ⋅ q1
    A2 = q2 ⋅ q1
    B1 = A2       # q1 ⋅ q2
    B2 = q2 ⋅ q2
    C1 = q1 ⋅ q3
    C2 = q2 ⋅ q3
    denom = A1 * B2 - A2 * B1
    if denom > 0
        a = (C1 * B2 - C2 * B1) / denom
        b = (C2 * A1 - C1 * A2) / denom
        sol.xpop .= @. (1 - a - b) * x1 + a * x2 + b * x3
    end
    sol
end


function warn_fat_lines(τl)
    nfat = count(>(1e5), τl)
    if nfat > 0
        @warn "Extremely optically thick lines found: $nfat"
    end
    nfat
end
warn_fat_lines(sol::Solution) = warn_fat_lines(sol.τl)


function isconverged(niter::Integer, nthick::Integer, τsum::Real,
        miniter::Integer, tolerance::Real)
    conv_min = niter >= miniter
    conv_sum = nthick == 0 || τsum / nthick < tolerance
    conv_min && conv_sum
end
function isconverged(sol::Solution, rdf::RunDef)
    isconverged(sol.niter, sol.nthick, sol.τsum, sol.iter.min, rdf.tolerance)
end


function first_iteration!(sol::Solution, rdf::RunDef)
    # NOTE The old level populations do not exist on first iteration, so store
    # results after first calculation.
    init_radiative!(sol, rdf)
    step_collision!(sol, rdf)
    solve_rate_eq!(sol)
    floor_population!(sol, rdf)
    store_population!(sol)
    init_tau_tex!(sol, rdf)
    sol
end


function iterate_solution!(sol::Solution, rdf::RunDef)
    clear_rates!(sol)
    step_radiative!(sol, rdf)
    step_collision!(sol, rdf)
    store_population!(sol)
    if rdf.reduced
        solve_rate_eq_reduced!(sol, rdf)
    else
        solve_rate_eq!(sol)
    end
    floor_population!(sol, rdf)
    step_tau_tex!(sol, rdf)
    # check conditions to accelerate
    after_start  = sol.niter >= sol.iter.start
    none_pending = (sol.niter - sol.iter.start) % sol.iter.pending == 0
    if after_start && none_pending
        ng_accelerate!(sol)
    end
    sol
end


function solve!(sol::Solution, rdf::RunDef; reuse_result=false)
    validate(sol, rdf)
    if reuse_result
        reset!(sol, keep_result=true)
        sol.niter = 1
    else
        sol.niter = 1
        first_iteration!(sol, rdf)
    end
    converged = false
    while !converged
        if sol.niter >= sol.iter.max
            @warn "Did not converge in $(sol.iter.max) iterations."
            break
        end
        sol.niter += 1
        iterate_solution!(sol, rdf)
        converged = isconverged(sol, rdf)
        sol.niter == 2 && warn_fat_lines(sol)
    end
    @debug "τl"   sol.τl[1:6]
    @debug "tex"  sol.tex[1:6]
    @debug "xpop" sol.xpop[1:6]
    converged, sol
end

function solve(rdf::RunDef; iter=nothing)
    iter = iter isa Nothing ? IterationParams() : iter
    sol  = Solution(rdf, iter)
    converged, _ = solve!(sol, rdf)
    converged, sol
end


@doc raw"""
    calc_radiation_temperature(xnu, tex, τl, intens_bg)

Calculate the radiation peak temperature or the background-subtracted
line intensity in units of the equivalent radiation temperature in the
Rayleigh-Jeans limit. This is equivalent to an observed "main beam" antenna
temperature ``T_\mathrm{mb}`` (antenna temperature corrected for the aperture
efficiency, spill-over, and atmospheric opacity) if the emission completely
fills the beam (i.e., a beam filling fraction of 1).

```math
T_R = \frac{c^2}{2 k \nu^2} \left( I^\mathrm{em}_\nu - I^\mathrm{bg}_\nu \right)
```
"""
function calc_radiation_temperature(
        xnu::F, tex::F, τl::F, intens_bg::F) where F <: Real
    # Black-body intensity in temperature units at the given excitation
    # temperature.
    bnutex = THC * xnu^3 * (inv ∘ expm1)(FK * xnu / tex)
    # Emergent intensity
    ftau = exp(-τl)
    intens_em = intens_bg * ftau + bnutex * (F(1) - ftau)
    # Calculate radiation temperature
    t_rad = FK / (THC * xnu * xnu) * (intens_em - intens_bg)
    t_rad
end


calc_wave(freq) = CLIGHT / freq / 1e5  # um
calc_kkms(Δv_cgs, t_rad) = GAUSS_AREA * Δv_cgs * t_rad / 1e5
calc_ergs(Δv_cgs, t_rad, xnu) = FGAUSS * KBOLTZ * Δv_cgs * t_rad * xnu^3


function get_results(sol::Solution, rdf::RunDef; freq_min=-Inf,
        freq_max=Inf)
    @assert freq_min < freq_max
    l = rdf.mol.levels
    t = rdf.mol.transitions
    # Compute limits on output from transition frequencies.
    # FIXME Add more helpful error messages.
    @assert maximum(t.spfreq) > freq_min
    @assert minimum(t.spfreq) < freq_max
    indices = findall(x -> freq_min < x < freq_max, t.spfreq)
    @assert length(indices) > 0
    freq = t.spfreq[indices]
    # Views would need to be copied to create the DataFrame, so just allocate
    # here instead with a slice and use `copycols=false` in the DataFrame
    # constructor.
    τl   = sol.τl[indices]
    tex  = sol.tex[indices]
    xnu  = t.xnu[indices]
    iupp = t.iupp[indices]
    ilow = t.ilow[indices]
    eup  = t.eup[indices]
    q_up = l.qnum[iupp]
    q_lo = l.qnum[ilow]
    backi = rdf.bg.backi[indices]
    # Population levels, note that iupp and ilow correspond to the indices in
    # the original (un-sliced) array.
    x_up = sol.xpop[iupp]
    x_lo = sol.xpop[ilow]
    # Calculate result values
    wave = calc_wave.(freq)
    t_rad = calc_radiation_temperature.(xnu, tex, τl, backi)
    flux_kkms = calc_kkms.(rdf.deltav, t_rad)
    flux_ergs = calc_ergs.(rdf.deltav, t_rad, xnu)
    # Convert density dictionary into a NamedTuple for splatting into DataFrame
    # constructor as separate columns.
    density = (; (Symbol("n_"*k) => v for (k,v) in rdf.density)...)
    DataFrame(
        trans_id=indices,
        q_up=q_up,
        q_lo=q_lo,
        e_up=eup,
        freq=freq,
        wave=wave,
        t_ex=tex,
        tau=τl,
        t_rad=t_rad,
        pop_up=x_up,
        pop_lo=x_lo,
        flux_kkms=flux_kkms,
        flux_ergs=flux_ergs,
        # Input parameters
        t_kin=rdf.tkin,
        ncolumn=rdf.cdmol,
        deltav=rdf.deltav/1e5,
        geometry=String(Symbol(rdf.escprob));
        density...,
        copycols=false,
     )
end
get_results(rdf::RunDef; kwargs...) = get_results(solve(rdf)[2], rdf; kwargs...)


function run(sol::Solution, rdfs; min_freq=-Inf, max_freq=Inf)
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


end  # module
