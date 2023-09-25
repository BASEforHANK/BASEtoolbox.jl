# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module BfH_EconFunc
# Sibling modules
using ..BfH_Tools

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
include("BfH_EconFunc/fcn_incomes.jl")
include("BfH_EconFunc/DistrSummaries.jl")
include("BfH_EconFunc/fcn_util_etc.jl")
end # module BASEforHANK
