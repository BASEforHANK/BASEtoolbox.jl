# General structure of examples

## Example folder structure

The folder of each example typically consists of the following:

`Model/`: Contains the model files specified by the user. These include:
- `input_aggregate_model.mod`: contains the user's model equations.
- `input_aggregate_names.jl`: contains the naming of the variables, in particular, which variables are states, choices, and shocks is specified here.
- `input_aggregate_steady_state.mod`: contains the equations for the steady state of the model.
- `input_compute_args_hh_prob_ss.jl`: provides the function that computes the arguments of the household problem in the steady state.
- `input_functions.jl`: contains the functions that can be user specified. CAN THE INPUTS BE CHANGED?
- `input_parameters.jl`: specifies the parameters. Keep in mind that you cannot take out parameters for the household model part.

`Data/`: Contains some data used in the baseline for the estimation. If you are not interested in estimation, this folder is not strictly required.

`main.jl`: The main file that runs the example and calls the methods of the toolbox. Occasionally this file is called `main_noestim.jl` to explicitly indicate that the corresponding model is not estimated.

## Setup of the model

The following guides through the files used to setup the aggregate model. They are located in the `Model/` folder of any example.

### Setup of the variables and names

The first step in setting up the model is to define the variables and their names. This is done in the `input_aggregate_names.jl` file. The user should usually specify or simply adjust only the following objects:

- `shock_names`: list of aggregate shocks as symbols (without any duplication).
- `state_names`: list of aggregate states as strings (without duplication). These names should be the same as the variable names used in `input_aggregate_model.mod`.
- `control_names`: list of aggregate controls as strings (without duplication). These names should be the same as the variable names used in `input_aggregate_model.mod`.
- `args_hh_prob_names`: list of arguments of the household problem as strings. These names should be the same as the variable names used in `input_aggregate_model.mod`. They are a subset of the variables defined above.
    !!! danger
    Importantly, some variables must exist in this list and with exactly these names as the household problem cannot be solved without them being specified. Naturally, those variables correspond to the variables showing up in the household's budget constraint like prices and taxes. See also the in-detail explanation of the [household problem](HouseholdProblem.md).

!!! tip
    The following objects need to be changed only if the multisector/-country feature is used:

- `n_rep`: integer specifying the duplicates, e.g. number of sectors or countries.
- `dup_state_names`: list of aggregate states as strings that should be duplicated.
- `dup_control_names`: list of aggregate controls as strings that should be duplicated.

!!! danger
    The following objects need to be changed only if changes outside the example structure are made:

- `distr_names`: list of cross-sectional controls or distributional summary variables as strings, which are also controls (but not listed there).

!!! note
    These are only the aggregate variables so that the distribution (state) and the marginal utilities/value functions (controls) are not included here.

### Setup of the steady state

The setup of the steady state distinguishes between two kinds of variables: household problem arguments and aggregate variables. The former are the variables the household needs to know to solve her individual planning problem, usually the variables that enter her budget constraint.

The set of variables comprising the household problem arguments are the ones specified in `args_hh_prob_names` in `input_aggregate_names.jl`. They are computed separately from the other steady state values (see below) in `input_compute_args_hh_prob_ss.jl`. They must be computable based on the capital stock in steady state and the model parameters. (Therefore those are the inputs of the function [`BASEforHANK.IncomesETC.compute_args_hh_prob_ss()`](@ref).) To facilitate computations, the user can use functions provided in `input_functions.jl`.

The steady state of all other aggregate variables is defined in `input_aggregate_steady_state.mod`. The syntax here is to refer to the variables as specified in `input_aggregate_names.jl` and add `$SS` to the variable name to refer to its steady state value. Once defined, a variable's steady state can be used for further definitions as if it is a value. Also here, the user can make use of functions provided in `input_functions.jl`.

!!! tip
    The command `@R$1` need to be changed only if the multisector/-country feature is used. In this case, the number should be set to the number of sectors/countries.

```@docs
BASEforHANK.IncomesETC.compute_args_hh_prob_ss
```

### More details on `input_functions.jl`

The file `input_functions.jl` includes multiple functions. Some of these are optional whereas others are used as part of the package. The functions are
 - `output()`: Calculates firm's total output from production inputs.
 - `profits_E_ss()`: Uses total output of the model economy to calculate the entrepreneurs' profits in steady state, which enter the budget considerations of the households.
 - `interest()`: Defines the net real interest rates, which households receive, also net of depreciation. In the baseline model it is equal to the marginal product of capital minus depreciation.
 - `wage()`: These are real wages, which are paid by the firm. Importantly these are not wages received by the household, which connects to the union profits, specified in the function below.
 - `profits_U_()`: Function, which use the difference in the wages paid by the firm and the wages paid to households to calculate union profits.

Please note that using these functions at the right places in your model equations in `input_aggregate_model.mod` ensures consistency of your model equations with the functions' specifications.

!!! warning
    Due to some current functionality, the input arguments to the functions `output()`, `wage()`, `interest()`, `profits_E_ss()` and `profits_U_()` cannot be changed. Keep this in mind when changing the aggregate model.

### Setup of the model parameters

The parameters are defined in the file `input_parameters.jl`. We discuss the structure concerning the parameters in more detail in the section [Parameters](#Structuring-parameters).

### Setup of the aggregate model equations

Finally, the file `input_aggregate_model.mod` contains the aggregate model equations. The user has to provide the equations in the form of
```julia
F[equation_number] = (lhs) - (rhs)
```
where the `equation_number` is based on a variable-specific index, following the pattern `indexes.<variable_name>`. The `<variable_name>` must be unique, i.e. each `equation_number` is used only once.

!!! danger
    The variables defined in `args_hh_prob_names` must exist as variables throughout the model, i.e. they must be defined (though possibly implicitly) in `input_aggregate_model.mod`.

#### Lags and leads

The toolbox comes with support only for future variables ("Prime"), i.e. when solving the model, they are automatically defined. To use them when writing down the equations, the user adds the suffix `Prime` to the variable name. So for example, the equation for capital accumulation, $K_{t+1} = (1-\delta_0)K_t + I_t$, is written as:
```julia
# Capital accumulation equation
F[indexes.I] = (log(KPrime)) - (log(K * (1.0 - m_par.δ_0) + I))
```

Of course, the user can also use lags by previously defining them accordingly. For example, to define real wage inflation, $\frac{w^F_t}{w^F_{t-1}} = \frac{\pi^w_t}{\pi_t}$, the user would define a variable `wFlag` first by using `wFlagPrime` and then use this new variable in the equation itself:
```julia
# Definition of lagged real wage
F[indexes.wFlag] = (log(wFlagPrime)) - (log(wF))

# Definition of real wage inflation
F[indexes.πw] = (log(wF / wFlag)) - (log(πw / π))
```

#### Auxiliary variables

You can specify auxiliary variables inside `input_aggregate_model.mod` by regular Julia code. They go without an equation number:

```julia
# Auxiliary variable

# Mass of households in each productivity state, distribution is (nb, nk, nh)
distr_h = sum(distrSS, dims = (1,2))
```

#### Setting variables entering the household problem

As mentioned above, the variables that enter the household problem must be defined in `input_aggregate_model.mod`. If the user does not want to use them, they can be set to a constant:
```julia
# Constant idiosyncratic income risk
F[indexes.σ] = (log(σ)) - (XSS[indexes.σSS])
```

#### Closing the aggregate model

In some (aggregate) equations, the households' decisions enter, e.g. market clearing conditions. They must be stated also in `input_aggregate_model.mod`, but are replaced in [`BASEforHANK.PerturbationSolution.Fsys()`](@ref) by the right aggregation during the preprocessing. However, for the estimation, they are used and are indeed correct as for the estimation only the derivative w.r.t. aggregates is important.

#### Special case: Multiple economies or sectors

The structure of the problem at hand allows for multiple economies or sectors in the aggregate model part as it would be nothing else than stating additional equations for each of them. To simplify writing up such models, the preprocessor allows using a certain syntax to mark the equations that would need to be repeated or copied for each economy and numbered accordingly. Instead, the preprocessor _automatically_ copies marked equations. Because many ASCII characters are reserved for other purposes, the syntax uses the `?` symbol.

The HANK economy is considered to be economy 1, not receiving an additional number. Any repeated equation will then get an additional number, starting with 2. The preprocessor will then automatically replace the `?` symbol with the respective number and will repeat the equation accordingly.

E.g., the following two code snippets from `input_aggregate_model.mod` are equivalent for the preprocessor:

```julia
# TFP HANK economy
F[indexes.Z]    = (log(ZPrime)) - (m_par.ρ_Z * log(Z))

# TFP economy 2
F[indexes.Z2]   = (log(Z2Prime)) - (m_par.ρ_Z * log(Z2))
```
```julia
# TFP HANK and all other economies
F[indexes.Z?]    = (log(Z?Prime)) - (m_par.ρ_Z * log(Z?))
```

At the beginning of the `input_aggregate_model.mod` file, the user has to specify the number of economies or sectors in their model by the following line:
```julia
# Setting the number of economies to 2.
@R?2
```

The same logic applies when specifying the variable names in `input_aggregate_names.jl` and `input_aggregate_steady_state.mod`. Consequently, the syntax can be also used when referring to steady state variables during writing down the model. E.g., ```Z?SS``` refers to the steady state value of variable `Z` for each economy.

In order to enable convenient syntax highlighting and because the `?` character serves other purposes in `Julia`, `.mod` files are used instead of `.jl` files for the user-specified model equations. In VS Code, you can still set the syntax highlighting to Julia for `.mod` files by clicking on the language mode in the bottom right corner of the editor.


## Structuring parameters

The file `input_parameters.jl` contains three structures to provide 1.) the model parameters, 2.) the numerical parameters, and 3.) the estimation settings.

### Model parameters

Each element of the `struct` [`ModelParameters`](@ref) consists of a parameter name, its value, its ascii name, its long name, its LaTeX name, its prior (please see the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl)-package for available options), and a Boolean whether the parameter should be estimated or not. If the parameter is not estimated, a prior may be suppressed by typing `_` instead of the prior.

Model parameters that only affect the aggregate dynamics can be freely adjusted.

!!! danger
    Parameters that (also) affect the household problem **have to** exist (with exactly these names and definitions). The reason for that is that these variables play a crucial role in the household problem, see also the in-detail explanation of the [household problem](HouseholdProblem.md).


### Numerical parameters

The `struct` [`NumericalParameters`](@ref) contains several types of parameters such as grids and meshes, variable types, and other parameters that determine the numerical approximation or solution technique, like `reduc` or `sol_algo`.

The user can modify some of the parameters and default values in `input_parameters.jl`. They serve as basis for [`call_find_steadystate()`](@ref), which sets other parameters automatically and adds further parameters from the discretization results like the transition matrix for productivity and the joint distribution.

In particular, `nh`, `nk`, and `nb` control the resolution for the individual productivity, illiquid asset, and liquid asset grid. The resolution of the copula used in the linearization does not need to coincide with that grid and is controlled by `nh_copula`, `nk_copula`, and `nb_copula`, respectively. Note, however, that the copula resolution should not exceed the actual grid size.


### Estimation settings

The `struct` [`EstimationSettings`](@ref) contains the settings for the estimation.

See also the section [Settings](Estimation.md#Settings).

## The main file

### Header

!!! tip
    Usually, the user does not need to change the header of the main file. The location of the files of the toolbox and the model are usually inferred from the path of the (activated) environment. See also [Installation](#index.md)

In this section, we set up the paths and pre-process the user's model inputs for the current example in the block
```julia
include(paths["src"] * "/Preprocessor/PreprocessInputs.jl");
```

Notably, this already pre-processes the aggregate model equations as well as the steady state file that are located in the folder `examples/<your_example>/Model/` and produces the generated functions in the folder `bld/<your_example>/Preprocessor/generated_fcns/`. For more details on the pre-processing, see [...].

Importantly, this pre-processing has to be performed before loading the BASEforHANK module defined in `BASEforHANK.jl`, which is then loaded via
```julia
include("BASEforHANK.jl")
using .BASEforHANK
```

`BASEforHANK.jl` is the key module file as it loads in the code base, sets up structures, and exports a number of functions and macros.

### Parameter and estimation preparation

The next step is to set up the model parameters, done by
```julia
m_par = ModelParameters();
e_set = BASEforHANK.e_set;
```

### Using the toolbox methods

After the steps described above, you are ready to use the toolbox' methods. A typical example could look as follows:

```julia
# Compute the steady state
ss_full = call_find_steadystate(m_par);

# Prepare the linearization, including the sparse DCT representation
sr_full = call_prepare_linearization(ss_full, m_par);

# Linearize the full model, i.e. find sparse state-space representation
lr_full = linearize_full_model(sr_full, m_par);
```

For more details, consider the complete example main files or the overview of [Methods](#index.md) provided by the toolbox.

## List of examples

We implemented the following examples. They may be a good starting point for your own model.

- The [baseline example](examples/baseline.md). This is a two-asset model. Additionally, we provide a one-asset and a complete-markets version of the baseline example.
- A [simpler version](examples/simpler_hank.md) of the baseline example: two assets, but simpler aggregate model part.
- A Krusell-Smith like model: one asset model with relatively simple aggregate model part.
- A [two-sector model](examples/two_sectors.md) with a service sector with monopolistic competition and nominal rigidities and a housing sector with a competitive producer.
