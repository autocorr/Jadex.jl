using Profile
using InteractiveUtils

using BenchmarkTools

using Jadex
using Jadex.Solver


function get_test_data(;reduced=false)
    datadir = joinpath(@__DIR__, "data")
    mol = Specie("hco+", datadir=datadir)
    bg  = blackbody_background(mol, tbg=2.730)
    #bg  = galactic_isrf(mol)
    rdf = RunDef(mol, density=Dict("h2" => 1e4), tkin=20.0, cdmol=1e13,
                 deltav=1.0, escprob=βsphere, bg=bg, reduced=reduced)
    sol = Solution(rdf)
    mol, bg, rdf, sol
end


function get_before_first_ldiv()
    _, _, rdf, sol = get_test_data()
    Solver.init_radiative!(sol, rdf)
    Solver.step_collision!(sol, rdf)
    sol
end


function get_first_iter()
    _, _, rdf, sol = get_test_data()
    Solver.first_iteration!(sol, rdf)
    sol
end


function get_solved(reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    niter, converged = solve!(sol, rdf)
    @info "niter=$niter, converged=$converged"
    sol
end


function bench_solution(reduced=false)
    _, _, rdf, _ = get_test_data(reduced=reduced)
    @benchmark begin
        sol = Solution(rdf)
        solve!(sol, rdf)
        df = get_results(sol, rdf; freq_max=500)
    end setup=(rdf=$rdf)
end


function profile_solution(reduced=false)
    mol, _, rdf, sol = get_test_data(reduced=reduced)
    @profile begin
        for _ in 1:1000
            solve!(sol, rdf)
            get_results(sol, rdf)
            reset!(sol)
        end
    end
    Profile.print()
end


function time_steps()
    _, _, rdf, sol = get_test_data(reduced=false)
    println(":: Initial iteration.")
    print("- init_rad:  "); @time Solver.init_radiative!(sol, rdf)
    print("- step_col:  "); @time Solver.step_collision!(sol, rdf)
    print("- solve:     "); @time Solver.solve_rate_eq!(sol)
    print("- floor_pop: "); @time Solver.floor_population!(sol, rdf)
    print("- store_pop: "); @time Solver.store_population!(sol)
    print("- init_tau:  "); @time Solver.init_tau_tex!(sol, rdf)
    print("- relax_pop: "); @time Solver.relax_population!(sol)
    println("\n:: Stepped iteration.")
    print("- clear:     "); @time Solver.clear_rates!(sol)
    print("- step_rad:  "); @time Solver.step_radiative!(sol, rdf)
    print("- step_col:  "); @time Solver.step_collision!(sol, rdf)
    print("- store_pop: "); @time Solver.store_population!(sol)
    rsol = deepcopy(sol)
    print("- solve_red: "); @time Solver.solve_rate_eq_reduced!(rsol, rdf)
    print("- solve:     "); @time Solver.solve_rate_eq!(sol)
    print("- floor_pop: "); @time Solver.floor_population!(sol, rdf)
    print("- step_tau:  "); @time Solver.step_tau_tex!(sol, rdf)
    print("- relax_pop: "); @time Solver.relax_population!(sol)
    println("\n:: Results.")
    print("- get_res:   "); @time Solver.get_results(sol, rdf)
    return nothing
end


function check_types()
    mol, _, rdf, sol = get_test_data()
    # initial iteration
    @code_warntype Solver.init_radiative!(sol, rdf)
    @code_warntype Solver.step_collision!(sol, rdf)
    @code_warntype Solver.solve_rate_eq!(sol)
    @code_warntype Solver.store_population!(sol)
    @code_warntype Solver.init_tau_tex!(sol, rdf)
    @code_warntype Solver.relax_population!(sol)
    # stepped iteration
    @code_warntype Solver.clear_rates!(sol)
    @code_warntype Solver.step_radiative!(sol, rdf)
    @code_warntype Solver.step_collision!(sol, rdf)
    @code_warntype Solver.store_population!(sol)
    @code_warntype Solver.solve_rate_eq!(sol)
    @code_warntype Solver.step_tau_tex!(sol, rdf)
    @code_warntype Solver.relax_population!(sol)
    return nothing
end

