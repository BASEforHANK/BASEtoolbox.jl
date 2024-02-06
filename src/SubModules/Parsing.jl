# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------


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
import FieldMetadata: label

export  ModelParameters,
        NumericalParameters,
        EstimationSettings,
        SteadyResults,
        LinearResults,
        EstimResults,
        SteadyStateStruct,
        IndexStruct,
        IndexStructAggr,
        produce_indexes,
        produce_indexes_aggr,
        metaflatten,
        prior,
        e_set,
        label

 export @writeXSS,
        @generate_equations,
        @make_fn,
        @make_fnaggr

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
