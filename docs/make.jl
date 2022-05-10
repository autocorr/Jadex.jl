using Documenter
include(joinpath(@__DIR__, "../src/Jadex.jl"))

makedocs(
    modules = [Jadex],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "Jadex",
    pages = [
             "index.md",
             "install.md",
             "userguide.md",
             "lib/api.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/autocorr/Jadex.jl.git",
    target = "build",
    forcepush = true,
)
