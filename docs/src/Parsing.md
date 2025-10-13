# Parsing
The `Parsing` submodule provides macros used for automatic code generation, structs, indexes, etc. These map, e.g., deviations from steady state of aggregate variables to equations lines in the non-linear difference equation that describe the economy.

## Collect variables: macros
```@docs
BASEforHANK.Parsing.@writeXSS
BASEforHANK.Parsing.@make_fn
BASEforHANK.Parsing.@make_fnaggr
BASEforHANK.Parsing.@make_struct
BASEforHANK.Parsing.@make_struct_aggr
BASEforHANK.Parsing.@generate_equations
BASEforHANK.Parsing.@write_args_hh_prob_ss
BASEforHANK.Parsing.@write_args_hh_prob
BASEforHANK.Parsing.@read_args_hh_prob
BASEforHANK.Parsing.@read_args_hh_prob_ss
```
