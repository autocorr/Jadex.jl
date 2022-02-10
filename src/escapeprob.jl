module EscapeProbability

export βsphere, βlvg, βslab


"""
    βsphere(τ)

Uniform sphere formula from Osterbrock (Astrophysics of Gaseous Nebulae and
Active Galactic Nuclei) Appendix 2 with power law approximations for large and
small optical depth.
"""
function βsphere(τ::Float64)
    τr = τ / 2.0
    if abs(τr) < 0.1
        return 1 - 0.75 * τr + τr^2 / 2.5 - τr^3 / 6 + τr^4 / 17.5
    elseif abs(τr) > 50
        return 0.75 / τr
    else
        a = 1 - 1 / 2τr^2
        b = 1 / τr + 1 / (2τr^2)
        return 0.75 / τr * (a + b * exp(-τ))
    end
end
βsphere(τ::Real) = βsphere(convert(Float64, τ))


"""
    βlvg(τ)

Expanding sphere, Large Velocity Gradient, or Sobolev case. Formula from De
Jong, Boland, and Dalgarno (1980, A&A 91, 68)[^DeJong80].  Corrected by a
factor of 2 in order to return 1 for τ=1.

[^DeJong80]: [ADS abstract](https://ui.adsabs.harvard.edu/abs/1980A%26A....91...68D/abstract)
"""
function βlvg(τ::Float64)
    τr = τ / 2.0
    if abs(τr) < 0.01
        return 1.0
    elseif abs(τr) < 7.0
        return 2.0 * (1.0 - exp(-2.34τr)) / 4.68τr
    else
        return 2.0 / (4.0τr * √log(τr / √π))
    end
end
βlvg(τ::Real) = βlvg(convert(Float64, τ))


"""
    βslab(τ)

Slab geometry (e.g. shocks) with power law approximations from De Jong,
Dalgarno, and Chu (1975, ApJ 199, 69)[^DeJong75].

[^DeJong75]: [ADS abstract](https://ui.adsabs.harvard.edu/abs/1975ApJ...199...69D/abstract)
"""
function βslab(τ::Float64)
    if abs(3τ) < 0.1
        return 1.0 - 1.5 * (τ + τ^2)
    elseif abs(3τ) > 50.0
        return 1.0 / 3τ
    else
        return (1.0 - exp(-3τ)) / 3τ
    end
end
βslab(τ::Real) = βslab(convert(Float64, τ))


end  # module
