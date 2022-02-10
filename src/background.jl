module Background
# TODO
# - Composition for adding fields together
# - OH dust
# NOTE Constructing fields all require knowing the `Specie` to evaluate
# the values at the transition frequencies, but seems redundant to pass
# Specie for the backgrounds and to the Solution. But if passing functions,
# then would have to re-evaluate each time which would be inefficient.

export BackgroundField, blackbody_background, galactic_isrf

using ..Constants: FK, THC, TCMB, EPS
using ..ReadData: Specie


"""
    BackgroundField(trj, backi, totalb)

# Arguments
- `trj`: Background brightness temperature `[K]` in the Rayleigh-Jeans limit.
- `backi`: Background intensity (integrated, ``J_ν``) in `[erg / (s cm^2 Hz sr)]`.
- `totalb`: Total background intensity. Note that because internal radiation
            fields are not implemented, `totalb` should be identical to `backi`.
"""
struct BackgroundField{F <: AbstractFloat}
    trj::Vector{F}
    backi::Vector{F}
    totalb::Vector{F}
end


"""
    planck(ΔE, tex)

Planck function to compute black body intensity in `[erg / (s cm^2 Hz sr)]`.

# Arguments
- `ΔE`: Energy level difference in wavenumber `[1/cm]`.
- `tex`: Excitation temperature in `[K]`.
"""
@inline function planck(ΔE, tex)
    hnu = FK * ΔE / tex  # unitless
    if hnu >= 160.0
        return EPS
    else
        return THC * ΔE^3 / expm1(hnu)
    end
end
planck(mol::Specie, tex) = planck(mol.transitions.xnu, tex)


"""
    blackbody_background(xnu; tbg::Real=TCMB)

Set the background radiation field from a blackbody with a given background
temperature. The default temperature is set for the CMB.
"""
function blackbody_background(xnu; tbg=TCMB)
    dtype = eltype(xnu)
    ntran = length(xnu)
    trj    = Array{dtype}(undef, ntran)
    backi  = Array{dtype}(undef, ntran)
    totalb = Array{dtype}(undef, ntran)
    backi .= planck.(xnu, tbg)
    trj .= tbg
    totalb .= backi
    BackgroundField(trj, backi, totalb)
end
function blackbody_background(mol::Specie; tbg=TCMB)
    blackbody_background(mol.transitions.xnu, tbg=tbg)
end


function galactic_isrf(xnu; kwargs...)
    backi = calc_galactic_isrf.(xnu; kwargs...)
    trj = @. FK * xnu / log1p(THC * xnu^3 / backi)
    totalb = copy(backi)
    BackgroundField(trj, backi, totalb)
end
galactic_isrf(mol::Specie) = galactic_isrf(mol.transitions.xnu)


@doc raw"""
    calc_galactic_isrf(xnu; uv_scaling=1.0)

Galactic interstellar radiation field (ISRF) as parametrized by Hocuk et al.
(2017; "H17") in Appendix B based on Zucconi et al. (2001). Six weighted
black-body components are used corresponding to the values in H17 Table B1.
The UV component of the ISRF is based on Draine (1978) and added according
to the polynomial in H17 Equation B2.

The angularly-integrated intensity ``J_ν`` is returned with units of
`erg / (s cm^2 sr Hz)`.
"""
function calc_galactic_isrf(xnu; uv_scaling=1.0)
    hnu = THC * xnu / 2  # erg
    return (
            # CMB component.
            1.0e+00 * planck(xnu, TCMB) +
            # Weighted black body intensities.
            1.0e-14 * planck(xnu, 7.50e3) +
            1.0e-13 * planck(xnu, 4.00e3) +
            4.0e-13 * planck(xnu, 3.00e3) +
            3.4e-09 * planck(xnu, 2.50e2) +
            2.0e-04 * planck(xnu, 2.33e1) +
            # Scaled polynomial for the UV-component.
            # NOTE Polynomial expects `hν` in [erg].
            uv_scaling * (
                4.28e03 * hnu^2 - 3.47e14 * hnu^3 + 6.96e24 * hnu^4
            )
    )
end


end  # module
