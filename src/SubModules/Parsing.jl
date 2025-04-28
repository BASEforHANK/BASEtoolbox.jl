# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module Parsing
# Sibling modules
using ..Tools
using ..Types

# 3rd Party modules
using LinearAlgebra,
    Parameters, Setfield, Flatten, FieldMetadata, LaTeXStrings, Distributions, Optim, Printf

import Flatten: flattenable

export ModelParameters,
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
    label,
    state_names,
    state_names_ascii,
    control_names,
    control_names_ascii,
    distr_names,
    shock_names,
    aggr_names,
    aggr_names_ascii,
    find_field_with_value,
    args_hh_prob_names,
    OneAsset,
    TwoAsset,
    AbstractMacroModel,
    mapround

export @writeXSS,
    @generate_equations,
    @make_fn,
    @make_fnaggr,
    @write_args_hh_prob,
    @read_args_hh_prob,
    @write_args_hh_prob_ss,
    @read_args_hh_prob_ss

# Documentation mode: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    include("../../examples/baseline/Model/input_aggregate_names.jl")
else
    include(Main.paths["src_example"] * "/Model/input_aggregate_names.jl")
end
# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
# Documentation mode: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    include("../../examples/baseline/Model/input_parameters.jl")
else
    include(Main.paths["src_example"] * "/Model/input_parameters.jl")
end
include("Parsing/Structs.jl")
include("../Preprocessor/prior.jl")

e_set = EstimationSettings(; shock_names = shock_names)
@make_struct IndexStruct
@make_struct_aggr IndexStructAggr
include("Parsing/MacroUtils.jl")

@make_fn produce_indexes
@make_fnaggr produce_indexes_aggr

include("Parsing/util.jl")

end # module Parsing
