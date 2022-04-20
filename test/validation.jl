"""
Note that for full validation of all LAMDA files available as of April 2022
that the following files should be corrected for improper conventions:
  * "c+" and "c+@uv" must have the trailing collision partners with preceding
    "!" removed at the end of the file.
  * "oatom@lique" must have the number of colliders amended from 5 to 6.
"""

using Test
using PyCall

using Jadex

DATADIR = ENV["RADEX_DATAPATH"]
RADEX = pyimport("spectralradex.radex")


function run_validation_tests()
    all_species = replace.(readdir(DATADIR), ".dat"=>"")
    @testset "Jadex Validation Tests" begin
        @testset "LAMDA file parsing" begin
            @testset "$name" for name in all_species
                @test Specie(name, datadir=DATADIR) !== nothing
            end
        end
        @testset "Execution tests" begin
            @testset "$name" for name in all_species
                mol = Specie(name, datadir=DATADIR)
                density = Dict(mol.colliders[1].name => 1e3)
                rdf = RunDef(mol, density=density)
                @test get_results(rdf) !== nothing
            end
        end
    end
end
