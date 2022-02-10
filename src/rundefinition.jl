module RunDefinition

export RunDef

using ..Constants: FK
using ..ReadData: Specie
using ..EscapeProbability: βlvg
using ..Background: BackgroundField, blackbody_background


struct RunDef{F <: AbstractFloat, B <: Function}
    mol::Specie{F}          # molecule container
    density::Dict{String, F}  # collision partner number densities, cm^-3
    totdens::F                # total number density of all partners, cm^-3
    tkin::F                   # kinetic temperature, K
    cdmol::F                  # molecular column density, cm^-2
    deltav::F                 # FWHM line width, cm s^-1
    escprob::B                # escape probability geometry
    bg::BackgroundField{F}    # radiation field background
    crate::Matrix{F}          # product of rate coeff and density
    ctot::Vector{F}           # total rate coeff for each level
    reduced::Bool             # reduce solution to coll. coupled
    nreduced::Int             # reduced number of levels
    minpop::F                 # minimum level population
end

function RunDef(
        mol::Specie{F};
        density=Dict("h2" => 1e4),
        tkin=20,
        cdmol=1e13,
        deltav=1.0,
        escprob::Function=βlvg,
        bg::Union{Nothing, BackgroundField{F}}=nothing,
        reduced::Bool=false,
        minpop=1e-20,
    ) where {F <: AbstractFloat}
    partner_names = [c.name for c in mol.colliders]
    @assert length(density) > 0
    @assert length(density) == mol.npart
    @assert all(p in partner_names for (p, n) in density)
    totdens = sum(n for (p, n) in density)
    @assert totdens > 0
    @assert tkin > 0
    @assert cdmol > 0
    @assert deltav > 0
    @assert minpop > 0
    deltav *= 1e5  # km/s to cm/s
    if bg === nothing  # use CMB blackbody
        bg = blackbody_background(mol)
    end
    density = Dict(p => F(d) for (p, d) in density)
    crate, ctot = get_collision_rates(mol, density, F(tkin))
    nreduced = count(<(10tkin/FK), mol.levels.eterm)
    if reduced && nreduced == mol.levels.n
        reduced = false
        @warn "All levels are collisionally coupled."
    end
    RunDef(mol, density, F(totdens), F(tkin), F(cdmol), F(deltav), escprob, bg,
           crate, ctot, reduced, nreduced, F(minpop))
end


"""
    get_collision_rates(mol::Specie, densities::Dict, tkin::Real)

Interpolate the collision rates for the kinetic temperature and sum over all
colliders.

# Returns
- `crate::Matrix`: Collision rate matrix (product of density and rate coeff.).
- `ctot::Vector`: Total collision rate for all transitions to a level.
"""
function get_collision_rates(mol::Specie, densities::Dict, tkin::Real)
    nlev = mol.levels.n
    partner_names = [c.name for c in mol.colliders]
    @assert all(p in partner_names for (p, n) in densities)
    colliders = Dict(c.name => c for c in mol.colliders)
    dtype = eltype(colliders[partner_names[1]].rates)
    crate = zeros(dtype, nlev, nlev)
    for (name, density) in densities
        collider = colliders[name]
        rate_temps = collider.temp
        tmin, tmax = extrema(rate_temps)
        rates = collider.rates
        # Compute downward collision rates.
        colld = zero(crate)
        # Extrapolate rates to lowest or highest values if the kinetic
        # temperature is beyond the extreme values.
        level_iter = enumerate(zip(collider.lcl, collider.lcu))
        # Clip to rates at lowest temperature.
        if tkin <= tmin
            tkin < tmin && @warn "Kinetic temperature lower than min value \
                                  for rates, $tkin < $tmin, using nearest."
            for (i, (e_lo, e_hi)) in level_iter
                colld[e_hi, e_lo] = rates[i, 1]
            end
        # Clip to rates at highest temperature.
        elseif tkin >= tmax
            tkin > tmax && @warn "Kinetic temperature greater than max value \
                                  for rates, $tkin > $tmax, using nearest."
            for (i, (e_lo, e_hi)) in level_iter
                colld[e_hi, e_lo] = rates[i, end]
            end
        # Interpolate between adjacent rates at the given temperature.
        else
            k_lo = findlast(rate_temps .< tkin)
            k_hi = findfirst(rate_temps .>= tkin)
            t_lo = rate_temps[k_lo]
            t_hi = rate_temps[k_hi]
            @assert k_lo !== nothing
            @assert k_hi !== nothing
            slope = (tkin - t_lo) / (t_hi - t_lo)
            for (i, (e_lo, e_hi)) in level_iter
                colld[e_hi, e_lo] =
                    rates[i,k_lo] + slope * (rates[i,k_hi] - rates[i,k_lo])
            end
        end
        @. crate += density * colld
    end
    # Compute upward rates from detailed balance.
    eterm = mol.levels.eterm
    gstat = mol.levels.gstat
    @inbounds for i_hi = 1:nlev, i_lo = 1:nlev
        ΔE = eterm[i_hi] - eterm[i_lo]
        etr = FK * ΔE / tkin
        if ΔE > 0
            # Avoid the exponential call and approximate the rate as zero if
            # the exponent is > 160 because exp(-160) ~ 1e-70
            if etr < 160
                g_hi = gstat[i_hi]
                g_lo = gstat[i_lo]
                crate[i_lo, i_hi] = g_hi / g_lo * exp(-etr) * crate[i_hi, i_lo]
            else
                crate[i_lo, i_hi] = zero(eltype(crate))
            end
        end
    end
    # Compute total collision rates from all levels to a given level.
    ctot = dropdims(sum(crate, dims=2), dims=2)
    crate, ctot
end


end  # module
