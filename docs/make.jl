push!(LOAD_PATH,"../src/")
using Documenter, BASEforHANK
makedocs(sitename="Documentation for BASEforHANK module",
            pages=[
                "Home" => "index.md",
                "Steady state" => "steadystate.md",
                "Linearization" => "linearization.md",
                "Estimation" => "estimation.md"
            ], format = Documenter.HTML(prettyurls = false))

 deploydocs(
    repo = "github.com/BASEforHANK/BASEtoolbox.jl.git", versions = nothing,
    )
