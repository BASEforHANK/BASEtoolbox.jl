# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------
# Packages used: Plots Distributions BenchmarkTools JLD2 FileIO DataFrames ForwardDiff
# SparseArrays LinearAlgebra Random LaTeXStrings MatrixEquations Roots KrylovKit JSON 
# CodecZlib SpecialFunctions FFTW Parameters Setfield MCMCChains StatsPlots Optim CSV 
# OrderedCollections Flatten FieldMetadata MKL

module Parsing
# Sibling modules
using ..Tools

# 3rd Party modules
using   LinearAlgebra,
        Parameters,
        Setfield,
        Flatten,
        FieldMetadata,
        LaTeXStrings,
        Distributions,
        Optim

import Flatten: flattenable    

export  ModelParameters,
        NumericalParameters,
        EstimationSettings,
        SteadyResults,
        LinearResults,
        EstimResults,
        SteadyState,
        IndexStruct,
        IndexStructAggr,
        produce_indexes,
        produce_indexes_aggr,
        metaflatten,
        prior,
        e_set

 export @writeXSS,
        @generate_equations

include("../Model/input_aggregate_names.jl")
# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("../Model/Parameters.jl")
include("Parsing/Structs.jl")
include("../Preprocessor/prior.jl")

e_set = EstimationSettings(shock_names = shock_names)
@make_struct IndexStruct
@make_struct_aggr IndexStructAggr
include("Parsing/MacroUtils.jl")


@make_fn produce_indexes
@make_fnaggr produce_indexes_aggr

end # module Parsing
