using PrecompileTools

@setup_workload begin
    datadir = joinpath(pkgdir(Jadex), "test", "data")
    @compile_workload begin
        βlvg(0.5)
        βslab(0.5)
        βsphere(0.5)
        mol = Specie("hco+", datadir=datadir)
        rdf = RunDef(mol, density=Dict("h2" => 1e4))
        sol = Solution(rdf)
        solve!(sol, rdf)
        reset!(sol)
        solve(rdf)
        get_results(rdf)
        get_results(sol, rdf)
    end
end
