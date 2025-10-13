using Documenter, BASEforHANK

makedocs(;
    sitename = "Documentation for BASEforHANK module",
    pages = [
        "Home" => "index.md",
        "Household problem" => "HouseholdProblem.md",
        "Computational Notes" => "ComputationalNotes.md",
        "Example Structure" => "GeneralStructure.md",
        "Steady state" => "SteadyState.md",
        "Perturbation solution" => "PerturbationSolution.md",
        "Estimation" => "Estimation.md",
        "Post estimation" => "PostEstimation.md",
        "Tools" => "Tools.md",
        "Parser" => "Parsing.md",
    ],
    format = Documenter.HTML(;
        mathengine = Documenter.HTMLWriter.MathJax3(),
        prettyurls = true,
    ),
)

deploydocs(; repo = "github.com/hildebrandecon/BASEtoolbox.jl.git", versions = nothing)
