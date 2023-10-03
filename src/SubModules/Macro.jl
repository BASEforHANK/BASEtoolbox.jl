# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module Macro

# Sibling modules
using ..Tools
using ..Parsing
using ..EconFunc
using ..Steady

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
include("Macro/compute_reduction.jl")
include("Macro/FSYS.jl")
include("Macro/LinearSolution.jl")
include("Macro/LinearSolution_estim.jl")
include("Macro/SolveDiffEq.jl")
include("Macro/Shuffle.jl")
include("Macro/update_model.jl")
include("Macro/model_reduction.jl")

include("../Preprocessor/generated_fcns/prepare_linearization_generated.jl")
include("../Preprocessor/generated_fcns/FSYS_agg_generated.jl")

end # module Macro
