"""
Note that for full validation of all LAMDA files available as of April 2022
that the following files should be corrected for improper conventions:
  * "c+" and "c+@uv" must have the trailing collision partners with preceding
    "!" removed at the end of the file.
  * "oatom@lique" must have the number of colliders amended from 5 to 6.
"""

using Test
using Printf
using PyCall
using DataFrames

using Jadex
using Jadex.Constants: TCMB

radex = pyimport("spectralradex.radex")

const DATADIR = ENV["RADEX_DATAPATH"]
const RADEX_COLLIDERS = Dict(
        "h2"   => 1,
        "p-h2" => 2,
        "o-h2" => 3,
        "e-"   => 4,
        "h"    => 5,
        "he"   => 6,
        "h+"   => 7,
)
const ZERO_DENSITY = Dict(c => 0.0 for c in keys(RADEX_COLLIDERS))


function get_radex_results(name, tkin, density, cdmol, linewidth=1.0; fmin=0.0,
            fmax=1e9, geometry=1)
    @assert in(geometry, [1,2,3])
    molfile = joinpath(DATADIR, "$name.dat")
    params = Dict(
            "molfile" => molfile,
            "tkin" => tkin,
            "tbg" => TCMB,
            "cdmol" => cdmol,
            "linewidth" => linewidth,
            "fmin" => fmin,
            "fmax" => fmax,
            "geometry" => geometry,
    )
    params = merge(ZERO_DENSITY, params, density)
    pd_df = radex.run(params)
    @assert pd_df !== nothing "RADEX Failed, check RADEX error messages."
    df = DataFrame()
    for col in pd_df.columns
        df[!, col] = getproperty(pd_df, col).values
    end
    rename!(df, Dict(
            :"E_UP (K)" => "e_up",
            :"WAVEL (um)" => "wave",
            :"T_ex" => "t_ex",
            :"T_R (K)" => "t_rad",
            :"POP UP" => "pop_up",
            :"POP LOW" => "pop_lo",
            :"FLUX (K*km/s)" => "flux_kkms",
            :"FLUX (erg/cm2/s)" => "flux_ergs",
            :"Qup" => "q_up",
            :"Qlow" => "q_lo",
    ))
    df
end


function get_dataframes(name, tkin, dens, cdmol)
    mol = Specie(name, datadir=DATADIR)
    first_partner = mol.colliders[1].name
    density = Dict(first_partner => dens)
    # Get results from Jadex
    rdf = RunDef(mol, density=density, tkin=tkin, cdmol=cdmol,
            escprob=βsphere)
    j_df = get_results(rdf)
    # Get results from Radex
    r_df = get_radex_results(name, tkin, density, cdmol)
    @assert size(r_df)[1] == size(j_df)[1]
    j_df, r_df
end


function max_residuals(name, tkin, dens, cdmol)
    j_df, r_df = get_dataframes(name, tkin, dens, cdmol)
    mask = r_df[!, "pop_up"] .> 1e-14 .&& 1e-3 .< r_df[!, "tau"] .< 1e1
    no_valid = !any(mask)
    function compare(col)
        j = j_df[mask, col]
        r = r_df[mask, col]
        δ = @. abs((j - r) / r)
        maximum(δ)
    end
    if no_valid
        return 0.0, 0.0, 0.0
    else
        δtex = compare("t_ex")
        δtau = compare("tau")
        δt_r = compare("t_rad")
        return δtex, δtau, δt_r
    end
end


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

        @testset "SpectralRadex integrity" begin
            params = radex.get_default_parameters()
            @test params["tkin"] ≈ 30.0
            @test radex.run(params) !== nothing
            df = radex.run(params)
            tau_r = [2.23009201e-04, 7.35152712e-04, 1.11196740e-03]
            @test all(df.tau.values[1:3] .≈ tau_r)
        end

        @testset "SpectralRadex comparison" begin
            ε = 1e-3
            all_tkin  = [  20,   40,   60,   80,  100]  # K
            all_dens  = [ 1e3,  1e4,  1e5,  1e6,  1e7]  # cm^-3
            all_cdmol = [1e11, 1e12, 1e13, 1e14]        # cm^-2
            test_species = ["13co", "co", "hco+", "hcn", "hnc", "n2hp", "p-nh3"]
            @testset "$name" for name in test_species
                @testset "$tkin, $dens, $cdmol" for
                            tkin  in all_tkin,
                            dens  in all_dens,
                            cdmol in all_cdmol
                    δtex, δtau, δt_r = max_residuals(name, tkin, dens, cdmol)
                    @test δtex .< ε
                    @test δtau .< ε
                    @test δt_r .< ε
                end
            end
        end
    end
end

