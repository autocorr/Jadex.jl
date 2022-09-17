using Base.Threads
using Profile
using InteractiveUtils

using BenchmarkTools

using Jadex
using Jadex.Solver
using Jadex.RunDefinition: get_collision_rates


function get_test_data(;reduced=false)
    datadir = joinpath(@__DIR__, "data")
    mol = Specie("hco+", datadir=datadir)
    bg  = blackbody_background(mol, tbg=2.730)
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


function get_solved(;reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    niter, converged = solve!(sol, rdf)
    @info "niter=$niter, converged=$converged"
    sol
end


function bench_collision_rates()
    mol, _, rdf, _ = get_test_data()
    density = Dict("h2" => 1e4)
    tkin = 20.0
    @benchmark get_collision_rates($mol, $density, $tkin)
end


function bench_solve(;reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    @benchmark begin
        solve!($sol, $rdf)
        reset!($sol)
    end
end


function bench_full_call(;reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    @benchmark begin
        solve!($sol, $rdf)
        df = get_results($sol, $rdf; freq_max=500)
        reset!($sol)
    end
end


function profile_collision_rates()
    mol, _, rdf, _ = get_test_data()
    density = Dict("h2" => 1e4)
    tkin = 20.0
    @profile begin
        for _ in 1:100_000
            get_collision_rates(mol, density, tkin)
        end
    end
    Profile.print()
end


function profile_solve(;reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    @profile begin
        for _ in 1:10_000
            solve!(sol, rdf)
            reset!(sol)
        end
    end
    Profile.print()
end


function profile_full_call(;reduced=false)
    _, _, rdf, sol = get_test_data(reduced=reduced)
    @profile begin
        for _ in 1:10_000
            solve!(sol, rdf)
            get_results(sol, rdf)
            reset!(sol)
        end
    end
    Profile.print()
end


function profile_big_grid()
    mol, bg, rdf, sol = get_test_data()
    solve!(sol, rdf)
    N = 10
    params = (
              exp10.(range(3, 5, N)),       # density
              collect(range(10, 30, N)),    # kinetic temperature
              exp10.(range(12, 13.5, N)),   # column density
              collect(range(0.5, 2.0, N)),  # linewidth
              [1, 2],                       # transitions
    )
    println("Number of threads: $(nthreads())")
    @profile begin
        rungrid(mol, params...; escprob=rdf.escprob, bg=rdf.bg)
    end
    Profile.print()
    @time rungrid(mol, params...; escprob=rdf.escprob, bg=rdf.bg)
    return nothing
end


function time_steps()
    _, _, rdf, sol = get_test_data(reduced=false)
    println(":: Escape probabilities.")
    print("- βsphere:   "); @time βsphere(1.0)
    print("- βlvg:      "); @time βlvg(1.0)
    print("- βslab:     "); @time βslab(1.0)
    println("\n:: Initial iteration.")
    print("- init_rad:  "); @time Solver.init_radiative!(sol, rdf)
    print("- step_col:  "); @time Solver.step_collision!(sol, rdf)
    print("- solve:     "); @time Solver.solve_rate_eq!(sol)
    print("- floor_pop: "); @time Solver.floor_population!(sol, rdf)
    print("- store_pop: "); @time Solver.store_population!(sol)
    print("- init_tau:  "); @time Solver.init_tau_tex!(sol, rdf)
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
    print("- ng_accel:  "); @time Solver.ng_accelerate!(sol)
    println("\n:: Results.")
    print("- get_res:   "); @time Solver.get_results(sol, rdf)
    return nothing
end


function check_types()
    mol, _, rdf, sol = get_test_data()
    # collision rates
    @code_warntype get_collision_rates(mol, rdf.density, rdf.tkin)
    # escape probability
    @code_warntype βsphere(1.0)
    @code_warntype βlvg(1.0)
    @code_warntype βslab(1.0)
    # initial iteration
    @code_warntype Solver.init_radiative!(sol, rdf)
    @code_warntype Solver.step_collision!(sol, rdf)
    @code_warntype Solver.solve_rate_eq!(sol)
    @code_warntype Solver.store_population!(sol)
    @code_warntype Solver.init_tau_tex!(sol, rdf)
    # stepped iteration
    @code_warntype Solver.clear_rates!(sol)
    @code_warntype Solver.step_radiative!(sol, rdf)
    @code_warntype Solver.step_collision!(sol, rdf)
    @code_warntype Solver.store_population!(sol)
    @code_warntype Solver.solve_rate_eq!(sol)
    @code_warntype Solver.step_tau_tex!(sol, rdf)
    @code_warntype Solver.ng_accelerate!(sol)
    return nothing
end

