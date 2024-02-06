# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module Estimation

# Sibling modules
using ..Tools
using ..Parsing
using ..IncomesETC
using ..SteadyState
using ..PerturbationSolution

# 3rd Party modules
using   LinearAlgebra,
        Distributions,
        Roots,
        Flatten,
        OrderedCollections,
        Setfield,
        Optim,
        FileIO,
        DataFrames,
        CSV,
        FieldMetadata,
        Parameters,
        FiniteDiff
    
using   Parameters: @unpack
using   MatrixEquations: lyapd
using   ProximalOperators: prox!, IndPSD
import  Flatten: flattenable

export  mode_finding,
        likeli,
        nearest_spd,
        rwmh 
        


include("Estimation/likeli.jl")
include("Estimation/filter_smoother.jl")
include("Estimation/mcmc.jl")
include("Estimation/measurement_error.jl")
include("Estimation/mode_finding.jl")
include("Estimation/nearest_spd.jl")

end # end submodule Estimation