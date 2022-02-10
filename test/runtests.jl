using Test
using Logging

using Jadex.Solver
using Jadex.ReadData
using Jadex.Background
using Jadex.RunDefinition
using Jadex.EscapeProbability


const TEST_DATA_DIR = "data"

Logging.disable_logging(Logging.Info)

prettyclose(x, y) = isapprox(x, y, rtol=1e-5)


"""
Test run configured to match the `example.inp` input file distributed with
RADEX and using the test LAMDA rate file for HCO+.
"""
function get_test_data(;reduced=false)
    mol = Specie("hco+", datadir=TEST_DATA_DIR)
    bg  = blackbody_background(mol, tbg=2.730)
    rdf = RunDef(mol, density=Dict("h2" => 1e4), tkin=20.0, cdmol=1e13,
                 deltav=1.0, escprob=βsphere, bg=bg, reduced=reduced)
    sol = Solution(rdf)
    mol, bg, sol, rdf
end


@testset "Jadex Tests" begin
    @testset "EscapeProbability" begin
        # Compare values to those produced by RADEX
        # By default `isapprox` compares to a relative tolerance
        # of √eps ~ 1.5e-8
        @test βsphere(1e-2) ≈ 0.996_259_979_202_38096
        @test βsphere(1e0)  ≈ 0.707_276_647_028_65400
        @test βsphere(1e2)  ≈ 1.499_700e-2

        @test βlvg(1e-3) ≈ 1.0
        @test isapprox(βlvg(1e0), 0.589_429_974_966_32554, rtol=1e-7)
        @test βlvg(1e2)  ≈ 5.472_036_662_316_83460e-3

        @test βslab(1e-2) ≈ 0.984_850
        @test βslab(1e0)  ≈ 0.316_737_643_877_37868
        @test βslab(1e2)  ≈ 3.333_333_333_333_33350e-3
    end

    @testset "ReadData" begin
        mol = Specie("hco+", datadir=TEST_DATA_DIR)
        @test mol.name == "HCO+"
        @test mol.levels.n == 21
        @test length(mol.levels.eterm) == mol.levels.n
        @test length(mol.levels.gstat) == mol.levels.n
        @test mol.transitions.n == 20
        @test mol.npart == 1
        collider = mol.colliders[1]
        @test collider.name == "h2"
        @test collider.ntran == 210
        @test collider.ntemp ==  12
        @test size(collider.rates) == (210, 12)
    end

    @testset "Rate interpolation" begin
        # Test interpolation at temperatures above and below defined values.
        # For HCO+ these are 10 and 400 K.
        mol = Specie("hco+", datadir=TEST_DATA_DIR)
        get_ctot_lo = () -> RunDef(mol, tkin=5).ctot[1]
        get_ctot_hi = () -> RunDef(mol, tkin=500).ctot[1]
        r_ctot_lo = 3.892_524_842_346_4167E-006
        r_ctot_hi = 3.226_799_565_449_8865E-005
        @test (@test_logs (:warn,) get_ctot_lo()) ≈ r_ctot_lo
        @test (@test_logs (:warn,) get_ctot_hi()) ≈ r_ctot_hi
    end

    @testset "Solver" begin
        # Test first iteration directly against values from RADEX. This is
        # mainly useful for testing that the rate matrix is solved correctly.
        @testset "HCO+ niter=1" begin
            mol, _, sol, rdf = get_test_data()
            Solver.init_radiative!(sol, rdf)
            Solver.step_collision!(sol, rdf)
            Solver.solve_rate_eq!(sol)
            r_rhs = [0.52942118255795823,       0.42348006153406509,
                     4.4332484887005752E-002,   2.4618512651810886E-003,
                     2.6271278329462406E-004,   3.5633840037103005E-005,
                     5.2622986435319052E-006,   7.2330917128867639E-007,
                     8.1653650236570583E-008,   5.5442857564631790E-009,
                     3.0872115628249801E-010,   1.7010498594352012E-011,
                     9.3571570860704663E-013,   3.8613906710454091E-014,
                     1.3375489717899483E-015,   3.2397859006308905E-017,
                     7.4739944934256157E-019,  -2.9408446823174254E-021,
                    -1.4530699175473601E-020,  -1.2538866586949176E-020,
                    -1.0744492622345223E-020,   1.0000000031710769E-026]
            @test sol.xpop[1:10] ≈ r_rhs[1:10]
        end
        # Test solution directly against values from RADEX.
        @testset "HCO+ rates" begin
            # Collision rate coefficients from RADEX.
            r_crate = [
                0.0000000000000000      5.5705964399852323E-006 3.1572656943602856E-006
                2.3000000000000000E-006 0.0000000000000000      4.1279993483040822E-006
                1.1999999999999999E-006 3.7999999999999996E-006 0.0000000000000000
            ]
            r_ctot = [1.1313904767896258E-005, 9.4209433386618675E-006,
                      9.8048174151457141E-006, 1.0776502314889678E-005,
                      1.1324231700306057E-005, 1.2336414444212895E-005,
                      1.2767931167545717E-005, 1.2544221514903589E-005,
                      1.3525836837778393E-005, 1.3029106509945849E-005,
                      1.3507607098126058E-005, 1.2698702656878599E-005,
                      1.2731742407370728E-005, 1.2257949734508990E-005,
                      1.2280634207361936E-005, 1.2251001413107053E-005,
                      1.1284951335437319E-005, 1.1370030871832578E-005,
                      1.1442545495723622E-005, 1.2528900484503935E-005,
                      1.1260799999999999E-005]
            mol, bg, sol, rdf = get_test_data()
            @test prettyclose(rdf.crate[1:3,1:3], r_crate)
            @test prettyclose(rdf.ctot, r_ctot)
        end
        @testset "HCO+ solve" begin
            # Solution results from RADEX using updated FGAUSS
            r_tex  = [4.5051798034847579,       3.7688070695763267,
                      3.7243181660539388,       6.0331342214233112,
                      9.7019419756071201,       12.272014457812368,
                      13.946518720012186,       14.791394998847824,
                      13.900297234400927,       14.324777007383883]
            r_xpop = [0.42207907860113225,      0.48965572583242456,
                      8.4188769278363940E-002,  3.7497174335928683E-003,
                      2.8231561650859021E-004,  3.8014026404015651E-005,
                      5.5429432689188009E-006,  7.4647575356132005E-007,
                      8.3596139005065791E-008,  5.8510927708706473E-009]
            r_taul = [4.6864228061150204,       5.3003854675683666,
                      0.88566341363879786,      3.6522458712078952E-002,
                      2.5266417455037499E-003,  3.2910148692731140E-004,
                      4.7722142804693704E-005,  6.4962868215111186E-006,
                      7.5106624178649448E-007,  5.2938138085011507E-008]
            # Solve full system and compare results to RADEX literals.
            mol, bg, sol, rdf = get_test_data()
            niter, converged = solve!(sol, rdf)
            @test converged
            @test niter == 29
            slice = 1:length(r_tex)
            @test prettyclose(sol.tex[slice],  r_tex)
            @test prettyclose(sol.xpop[slice], r_xpop)
            @test prettyclose(sol.τl[slice],   r_taul)
        end
        @testset "HCO+ reduced" begin
            # Solution results from RADEX using reduced=true and patched FGAUSS
            r_tex  = [4.5051798035302424,       3.7688070696279143,
                      3.7243181661972127,       6.0331342234490375,
                      9.7019419990226083,       12.272014638252440,
                      13.946520068600222,       14.791406148762562,
                      13.900420003449824,       13.963891350582442]
            r_xpop = [0.42207907860767147,      0.48965572584470751,
                      8.4188769283093420E-002,  3.7497174343008540E-003,
                      2.8231561683089257E-004,  3.8014026649782931E-005,
                      5.5429434752922248E-006,  7.4647593640220123E-007,
                      8.3596305332561226E-008,  5.8512475907139018E-009]
            r_taul = [4.6864228061592792,       5.3003854676823705,
                      0.88566341368469126,      3.6522458716810410E-002,
                      2.5266417467229687E-003,  3.2910148762991718E-004,
                      4.7722143271694613E-005,  6.4962871699427997E-006,
                      7.5106650828449555E-007,  5.3148449509563594E-008]
            # Solve reduced level system and compare.
            mol, bg, sol, rdf = get_test_data(reduced=true)
            niter, converged = solve!(sol, rdf)
            @test converged
            @test niter == 29
            slice = 1:length(r_tex)
            @test prettyclose(sol.tex[slice],  r_tex)
            @test prettyclose(sol.xpop[slice], r_xpop)
            @test prettyclose(sol.τl[slice],   r_taul)
        end
        @testset "HCO+ results" begin
            # Solution results from RADEX using patched FGAUSS and GAUSS_AREA
            # for the first five levels (1-0 to 5-4).
            r_t_r =  [1.5569810039583865,
                      0.59274361520339935,
                      0.17892104987462643,
                      3.7027833677387996e-2,
                      6.6645760438094120e-3]
            r_kkms = [1.6573550112012627,
                      0.63095606081085065,
                      0.19045556616632217,
                      3.9414909715100213e-2,
                      7.0942217507199805e-3]
            r_ergs = [1.5142715128148389e-8,
                      4.6117193143380676e-8,
                      4.6979363178822558e-8,
                      2.3043927507445803e-8,
                      8.1000449885733543e-9]
            # Solve and compare results
            mol, bg, sol, rdf = get_test_data()
            niter, converged = solve!(sol, rdf)
            df = get_results(sol, rdf, freq_max=500)
            @test prettyclose(df.t_rad, r_t_r)
            @test prettyclose(df.flux_kkms, r_kkms)
            @test prettyclose(df.flux_ergs, r_ergs)
        end
    end
end

