module Constants

export CLIGHT, HPLANCK, KBOLTZ, AMU, FK, THC
export TCMB
export FGAUSS, MINITER, MAXITER, CCRIT, EPS, MINPOP


# Physical constants in CGS
const CLIGHT  = 2.99792458e+10  # speed of light     (cm/s)
const HPLANCK = 6.62609630e-27  # Planck constant    (erg/Hz)
const KBOLTZ  = 1.38065050e-16  # Boltzmann constant (erg/K)
const AMU     = 1.67262171e-24  # atomic mass unit   (g)
const FK      = HPLANCK * CLIGHT / KBOLTZ
const THC     = 2 * HPLANCK * CLIGHT

# CMB Blackbody temperature from Fixsen (2009)
const TCMB = 2.72548  # +/- 0.00057 K

# Mathematical constants
# NOTE RADEX uses a truncated FWHM value of 1.0645 which
#      results in a value for `FGAUSS` of 26.753803007
const FGAUSS  = √π / (2 * √log(2)) * 8π

const MINITER = 10     # minimum number of iterations
const MAXITER = 9999   # maximum number of iterations

const CCRIT   = 1e-6   # relative tolerance of solution
const EPS     = 1e-30  # round-off error
const MINPOP  = 1e-20  # minimum level population


end  # module
