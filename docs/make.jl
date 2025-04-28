using Documenter, BASEforHANK

makedocs(;
    sitename = "Documentation for BASEforHANK module",
    pages = [
        "Home" => "index.md",
        "Steady state" => "SteadyState.md",
        "Perturbation solution" => "PerturbationSolution.md",
        "Estimation" => "Estimation.md",
        "Post estimation" => "PostEstimation.md",
        "Tools" => "Tools.md",
        "Parser" => "Parsing.md",
        "Household problem" => "HouseholdProblem.md",
        "Baseline example" => "examples/baseline.md",
        "Computational Notes" => "Computational Notes.md",
    ],
    format = Documenter.HTML(;
        mathengine = Documenter.HTMLWriter.MathJax3(),
        prettyurls = true,
    ),
)

deploydocs(; repo = "github.com/BASEforHANK/BASEtoolbox.jl.git", versions = nothing)
