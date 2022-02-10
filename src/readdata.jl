module ReadData

export VALID_PARTNERS
export EnergyLevels, Transitions, CollisionPartner, Specie
export show

using CSV
using DataFrames


const VALID_PARTNERS = ["h2", "p-h2", "o-h2", "e", "h", "he", "h+"]


struct EnergyLevels{F <: AbstractFloat}
    n::Int                    # number of energy levels
    eterm::Vector{F}          # energy levels, in cm^-1
    gstat::Vector{F}          # statistical weights
    qnum::Vector{String}      # quantum numbers
end


struct Transitions{F <: AbstractFloat}
    n::Int                    # number of radiative transitions
    iupp::Vector{Int}         # upper state
    ilow::Vector{Int}         # lower state
    aeinst::Vector{F}         # Einstein A
    spfreq::Vector{F}         # spectral line frequencies
    eup::Vector{F}            # upper state energy, E_u / k
    xnu::Vector{F}            # energy level difference
end


struct CollisionPartner{F <: AbstractFloat}
    name::String              # collision partner reference
    ntran::Int                # number of collisional transitions
    ntemp::Int                # number of collisional temperatures
    temp::Vector{F}           # temperatures
    lcu::Vector{Int}          # upper state index for collision trans.
    lcl::Vector{Int}          # lower state index for collision trans.
    rates::Matrix{F}          # collision rates, cm^3 s^-1
end


struct Specie{F <: AbstractFloat}
    name::String
    amass::F                  # specie weight
    levels::EnergyLevels{F}
    transitions::Transitions{F}
    npart::Int                # number of collision partners
    colliders::Vector{CollisionPartner}
end

@doc raw"""
    Specie(name::AbstractString; datadir=datadir)

Parse the LAMDA collision rate file for an atomic or molecular specie from the
filename stem (i.e., without the ".dat" file extension). For example, the stem
for the `"hco+@xpol.dat"` file would be `"hco+@xpol"`.
"""
function Specie(name::AbstractString; datadir="data")
    path = joinpath(datadir, "$name.dat")
    lines = strip.(readlines(path))
    # Parse header
    mname = String(lines[2])
    amass = parse(Float64, lines[4])
    @assert amass > 0
    # Parse energies
    levels = parse_energies(lines)
    # Parse transitions
    transitions = parse_transitions(lines, levels)
    # Parse colliders
    colliders = parse_colliders(lines, levels, transitions)
    npart = length(colliders)
    Specie(mname, amass, levels, transitions, npart, colliders)
end


"""
    parse_table(lines, types)

Parse a white-space delimited, fixed-width table into a set of columns with
known data types.
"""
function parse_table(lines, types)
    text = join(lines, "\n")
    file = CSV.File(
            IOBuffer(text),
            delim=' ',
            ignorerepeated=true,
            header=false,
            types=types,
    )
    DataFrame(file)
end


"""
    parse_energies(lines)

Energy level terms begin on the 8th line and are of the form:

```
!LEVEL + ENERGIES(cm^-1) + WEIGHT + J
   1     0.000000000   1.0     0
```
"""
function parse_energies(lines)
    nlev = parse(Int, lines[6])
    @assert nlev >= 1
    types = [Int, Float64, Float64, String]
    work_lines = lines[8:7+nlev]
    df = parse_table(work_lines, types)
    eterm = df.Column2
    gstat = df.Column3
    qnum  = df.Column4
    EnergyLevels(nlev, eterm, gstat, qnum)
end


"""
    parse_transitions(lines, levels::EnergyLevels)

Transition terms are of the form:

```
!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
    1     2     1  4.251e-05          89.18839570     4.28
```
"""
function parse_transitions(lines, levels::EnergyLevels)
    ntran = parse(Int, lines[9+levels.n])
    @assert ntran >= 1
    i_start = levels.n + 11
    work_lines = lines[i_start:i_start+ntran-1]
    types = [Int, Int, Int, Float64, Float64, Float64]
    df = parse_table(work_lines, types)
    iupp = df.Column2
    ilow = df.Column3
    aeinst = df.Column4
    spfreq = df.Column5
    eup = df.Column6
    xnu = levels.eterm[iupp] .- levels.eterm[ilow]
    Transitions(ntran, iupp, ilow, aeinst, spfreq, eup, xnu)
end


"""
    parse_colliders(lines, levels, transitions)

Collider collision rate terms are of the form:

```
!NUMBER OF COLL PARTNERS
1
!COLLISIONS BETWEEN
1  HCO+ - H2 from Flower (1999)
!NUMBER OF COLL TRANS
210
!NUMBER OF COLL TEMPS
12
!COLL TEMPS
   10.0   20.0   30.0   50.0   70.0  100.0  150.0  200.0  250.0  300.0  350.0  400.0
!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)
    1     2     1  2.6e-10 2.3e-10 2.1e-10 2.0e-10 1.9e-10 1.8e-10 2.0e-10 2.2e-10 2.3e-10 2.5e-10 2.7e-10 2.8e-10
```
"""
function parse_colliders(lines, levels, transitions)
    npart = parse(Int, lines[12+levels.n+transitions.n])
    icol = 0
    colliders = Array{CollisionPartner}(undef, npart)
    # Instead of counting lines, read lines until reaching the comment
    # "!COLLISIONS BETWEEN", "!Collisions with", etc.
    for (i, line) in enumerate(lines)
        if occursin(r"^!COLLISIONS"i, line)
            icol += 1
            name = VALID_PARTNERS[parse(Int, lines[i+1][1])]
            ntran = parse(Int, lines[i+3])
            ntemp = parse(Int, lines[i+5])
            temp = parse.(Float64, split(lines[i+7]))
            @assert ntemp == length(temp)
            types = vcat([Int, Int, Int], repeat([Float64], ntemp))
            work_lines = lines[i+9:i+8+ntran]
            df = parse_table(work_lines, types)
            @assert ntran == size(df, 1)
            lcu = df.Column2
            lcl = df.Column3
            rates = Matrix(df[:,4:end])
            collider = CollisionPartner(name, ntran, ntemp, temp, lcu, lcl, rates)
            colliders[icol] = collider
        end
    end
    if icol != npart
        error("Could not parse all colliders.")
    end
    colliders
end


end  # module
