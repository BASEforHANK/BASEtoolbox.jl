# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module PerturbationSolution

# Sibling modules
using ..Tools
using ..Parsing
using ..IncomesETC
using ..SteadyState

# 3rd Party modules
using   LinearAlgebra,
        SparseArrays,
        BlockDiagonals,
        Distributions,
        Roots,
        ForwardDiff,
        Flatten,
        Setfield
    
using   Parameters: @unpack
using   MatrixEquations: lyapd

export  LinearSolution,
        LinearSolution_estim,
        compute_reduction,
        prepare_linearization,
        model_reduction,
        update_model

        

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("../Model/input_aggregate_names.jl") # this dependency is not ideal
include("PerturbationSolution/compute_reduction.jl")
include("PerturbationSolution/FSYS.jl")
include("PerturbationSolution/LinearSolution.jl")
include("PerturbationSolution/LinearSolution_estim.jl")
include("PerturbationSolution/SolveDiffEq.jl")
include("PerturbationSolution/Shuffle.jl")
include("PerturbationSolution/update_model.jl")
include("PerturbationSolution/model_reduction.jl")

include("../Preprocessor/generated_fcns/prepare_linearization_generated.jl")
include("../Preprocessor/generated_fcns/FSYS_agg_generated.jl")

end # module BASEforHANK.PerturbationSolution
