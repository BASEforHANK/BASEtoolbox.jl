# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module Steady
# Sibling modules
using ..Tools
using ..Parsing
using ..EconFunc

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
include("../Model/input_aggregate_names.jl") # this dependency is not ideal
include("Steady/EGM/EGM_policyupdate.jl")
include("Steady/EGM/updateV.jl")
include("Steady/find_steadystate.jl")
include("Steady/first_stage_reduction.jl")
include("Steady/IM_fcns/fcn_directtransition.jl")
include("Steady/IM_fcns/fcn_kdiff.jl")
include("Steady/IM_fcns/fcn_ksupply.jl")
include("Steady/IM_fcns/fcn_makeweights.jl")
include("Steady/IM_fcns/fcn_maketransition.jl")
include("Steady/IM_fcns/fcn_MultipleDirectTransitions.jl")
end # module Steady
