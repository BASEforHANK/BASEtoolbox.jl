# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module BfH_Steady
# Sibling modules
using ..BfH_Tools
using ..BfH_Parsing
using ..BfH_EconFunc

# 3rd Party modules
using   LinearAlgebra,
        SparseArrays,
        Distributions,
        Roots,
        ForwardDiff,
        Flatten
    
using Parameters: @unpack
using KrylovKit: eigsolve
using FFTW: dct, ifft

export  updateV!,
        updateV,
        EGM_policyupdate!,
        EGM_policyupdate,
        DirectTransition!,
        DirectTransition,
        Ksupply,
        find_steadystate,
        first_stage_reduction
        

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("../1_Model/input_aggregate_names.jl") # this dependency is not ideal
include("BfH_Steady/EGM/EGM_policyupdate.jl")
include("BfH_Steady/EGM/updateV.jl")
include("BfH_Steady/find_steadystate.jl")
include("BfH_Steady/first_stage_reduction.jl")
include("BfH_Steady/IM_fcns/fcn_directtransition.jl")
include("BfH_Steady/IM_fcns/fcn_kdiff.jl")
include("BfH_Steady/IM_fcns/fcn_ksupply.jl")
include("BfH_Steady/IM_fcns/fcn_makeweights.jl")
include("BfH_Steady/IM_fcns/fcn_maketransition.jl")
include("BfH_Steady/IM_fcns/fcn_MultipleDirectTransitions.jl")
end # module BfH_Steady
