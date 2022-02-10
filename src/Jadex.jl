module Jadex

export βsphere, βlvg, βslab
export Specie
export BackgroundField, blackbody_background, galactic_isrf
export RunDef
export Solution, reset!, solve!, get_results

include("constants.jl")
include("wrap_slatec.jl")
include("readdata.jl")
include("escapeprob.jl")
include("background.jl")
include("rundefinition.jl")
include("solver.jl")

using .Solver
using .ReadData
using .Background
using .RunDefinition
using .EscapeProbability

end
