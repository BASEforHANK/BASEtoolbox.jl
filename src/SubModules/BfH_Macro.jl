# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module BfH_Macro

# Sibling modules
using ..BfH_Tools
using ..BfH_Parsing
using ..BfH_EconFunc
using ..BfH_Steady

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
        prepare_linearization
        

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("../1_Model/input_aggregate_names.jl") # this dependency is not ideal
include("BfH_Macro/compute_reduction.jl")
include("BfH_Macro/FSYS.jl")
include("BfH_Macro/LinearSolution.jl")
include("BfH_Macro/LinearSolution_estim.jl")
include("BfH_Macro/SolveDiffEq.jl")
include("BfH_Macro/Shuffle.jl")

include("../Preprocessor/generated_fcns/prepare_linearization_generated.jl")
include("../Preprocessor/generated_fcns/FSYS_agg_generated.jl")

end # module BfH_Macro
