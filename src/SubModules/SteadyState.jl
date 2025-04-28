# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module SteadyState
# Sibling modules
using ..Tools
using ..Parsing
using ..IncomesETC
using ..Types

# 3rd Party modules
using LinearAlgebra, SparseArrays, Distributions, Roots, ForwardDiff, Flatten, Printf

using Parameters: @unpack
using KrylovKit: eigsolve
using FFTW: dct, ifft

export updateW!,
    updateW,
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
# Documentation mode: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    include("../../examples/baseline/Model/input_aggregate_names.jl") # this dependency is not ideal
else
    include(Main.paths["src_example"] * "/Model/input_aggregate_names.jl") # this dependency is not ideal
end
include("SteadyState/EGM/EGM_policyupdate.jl")
include("SteadyState/EGM/updateW.jl")
include("SteadyState/find_steadystate.jl")
include("SteadyState/first_stage_reduction.jl")
include("SteadyState/IM_fcns/fcn_directtransition.jl")
include("SteadyState/IM_fcns/fcn_kdiff.jl")
include("SteadyState/IM_fcns/fcn_ksupply.jl")
include("SteadyState/IM_fcns/fcn_makeweights.jl")
include("SteadyState/IM_fcns/fcn_maketransition.jl")
include("SteadyState/IM_fcns/fcn_MultipleDirectTransitions.jl")
end # module BASEforHANK.SteadyState
