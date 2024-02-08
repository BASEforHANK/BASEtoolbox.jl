#cd("docs")
using Documenter
push!(LOAD_PATH,"../src/")
 #include("../src/BASEforHANK.jl")
 using BASEforHANK
makedocs(sitename="Documentation for BASEforHANK module",
            pages=[
                "Home" => "index.md",
               "Steady state" => "SteadyState.md",
                "Perturbation solution" => "PerturbationSolution.md",
                "Estimation" => "Estimation.md",
                "Post estimation" => "PostEstimation.md",
                "Tools" => "Tools.md",
                "Parser" => "Parsing.md"
            ], format = Documenter.HTML(prettyurls = false))

 deploydocs(
    repo = "github.com/BASEforHANK/BASEtoolbox.jl.git", versions = nothing,
    )
