# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module EconFunc
# Sibling modules
using ..Tools

# 3rd Party modules
using   LinearAlgebra,
        ForwardDiff

export  distrSummaries,
        incomes!,
        incomes,
        util,
        mutil!,
        mutil,
        invmutil!,
        invmutil,
        employment,
        interest,
        wage,
        output,
        profitsSS_fnc,
        qΠSS_fnc,
        value_liquid

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("EconFunc/fcn_incomes.jl")
include("EconFunc/DistrSummaries.jl")
include("EconFunc/fcn_util_etc.jl")
end # module BASEforHANK
