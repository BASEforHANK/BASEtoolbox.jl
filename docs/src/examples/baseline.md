# Baseline Example

This baseline example comes closest to the original paper ([Bayer, Born, and Luetticke (2024, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20201875)). We provide two mainboards that accompany this example in the `examples/baseline/` folder:

1. `main.jl` which showcases the toolbox including the estimation part and
2. `main_noestim.jl` which showcases the toolbox without the estimation part.

In the following, we focus on `main.jl`.

Please refer to the main page of the documentation for general information on how to set up the toolbox.

The provided `main.jl` shows how a typical estimation proceeds in three main steps. First, we solve the steady state of the model. Then, the algorithm performs a two-step dimensionality reduction as described in the accompanying paper ([Bayer, Born, and Luetticke (2024, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20201875)). The second step of this reduction uses the prior information to obtain and approximate factor representation from an initial, not further reduced solution. Secondly, we compute the linearized dynamics of the reduced model around the steady state. Thirdly, we construct the likelihood of the model parameters given the data and use Bayesian methods to estimate them. More details on the three steps are provided below and in the respective sections of the documentation.

Next, we move through `main.jl` step by step.

## Setup of the model

To define the aggregate part of the baseline example's model, we included the aggregate model block in `examples/baseline/Model/`.

Firstly, we define the aggregate variables in `input_aggregate_names.jl`. In particular, we define a (symbolic) list of *shocks* as `shock_names`, a (string) list of *state variables* as `state_names`, and a (string) list of *control variables* as `control_names`. We also define a (string) list of *distributional statistics* that are just controls as `distr_names`.

Recall: These are only the aggregate variables so that the distribution (state) and the marginal utilities/value functions (controls) are not included here.

Importantly, some variables **have to** exist (with exactly these names and definitions). The reason for that is that these variables play a crucial role in the household problem, see also the in-detail explanation of the [`HouseholdProblem.md`](../HouseholdProblem.md).

The steady state of the aggregate variables is defined in `input_aggregate_steady_state.mod`. Notably, there are a few exceptions which are provided using a different syntax. The reason for this is that these variables affect the household problem. Therefore, they are computed in `compute_args_hh_prob_ss` based on the user-provided functions in `input_functions.jl` as well as assumptions of the household problem and some further assumptions. If the user provides these steady state values also in `input_aggregate_steady_state.mod`, these will be compared to the values computed in `compute_args_hh_prob_ss` and a warning will be issued if they differ.

The file `input_parameters.jl` contains three structures to provide model parameters, numerical parameters, and estimation settings.

Model parameters that only affect the aggregate dynamics can be freely adjusted. Parameters that (also) affect the household problem **have to** exist (with exactly these names and definitions). The reason for that is that these variables play a crucial role in the household problem, see also the in-detail explanation of the [`HouseholdProblem.md`](../HouseholdProblem.md). Each element of the `struct` [`ModelParameters`](@ref) consists of a parameter name, its value, its ascii name, its long name, its LaTeX name, its prior (please see the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)-package for available options), and a Boolean whether the parameter should be estimated or not.

The model parameters need to be set in your mainboard.

The `struct` [`NumericalParameters`](@ref) numerical parameters contain several types of parameters such as grids and meshes, variable types etc.

The numerical parameters are automatically set in `call_find_steadystate()`, based on their default values in `input_parameters.jl`, where the user can modify them.

Finally, the file `input_aggregate_model.jl` contains the aggregate model equations. The user has to provide the equations in the form of
```julia
F[equation number] = (lhs) - (rhs)
```
where the equation number is based on an variable-specific index, to be explained later.


## Header

In this section, we set up the paths and pre-process the user's model inputs for the current example in the block
```julia
include(paths["src"] * "/Preprocessor/PreprocessInputs.jl");
```

Notably, this already pre-processes the aggregate model equations as well as the steady state file that are located in the folder `examples/baseline/Model/` and produces the generated functions in the folder `bld/baseline/Preprocessor/generated_fcns/`. For more details on the pre-processing, see [...].

Importantly, this pre-processing has to be performed before loading the BASEforHANK module defined in `BASEforHANK.jl`, that is then loaded via
```julia
include("BASEforHANK.jl")
using .BASEforHANK
```

`BASEforHANK.jl` is the key module file as it loads in the code base, sets up structures, and exports a number of functions and macros.


## Steady state and first dimensionality reduction
The command
```
sr_full = compute_steadystate(m_par)
```
calls the functions [`BASEforHANK.find_steadystate()`](@ref) and [`BASEforHANK.prepare_linearization()`](@ref) and saves their returns in an instance `sr_full` of the `struct` `SteadyResults`.
In exact, `sr_full` contains vectors of the steady-state variables (together with index-vectors to reference them by name) and
the steady-state distribution of income and assets. It also contains the marginal value functions and the distributions as well as their first-stage model reduction counterparts (obtained through DCTs).

!!! tip
    `sr_full` may be saved to the local file system by calling
    ```
    @save "Output/Saves/steadystate.jld2" sr_full
    ```
    and can be loaded for a future session with
    ```
    @load "Output/Saves/steadystate.jld2" sr_full
    ```
More details can be found in the section "Steady State".

## Linearize full model
After computing the steady state and saving it in the `SteadyResults`-struct named `sr_full`,
```
lr_full = linearize_full_model(sr_full, m_par)
```
computes the linear dynamics of the "full" model, i.e., using the first-stage model reduction, around the steady state (in the background, this calls [`BASEforHANK.PerturbationSolution.LinearSolution()`](@ref)) and saves a state-space representation in the instance `lr_full` of the `struct` `LinearResults` (see [`linearize_full_model()`](@ref)).

Linearization of the full model takes a few seconds. The resulting state space is relatively large, because the copula and the value functions are treated fully flexible in this first step. As a result, also computing the first-order dynamics of this model takes a few seconds as well.

## Model reduction
This large state-space representation can, however, be reduced substantially using an approximate factor representation. For this purpose, we run
```
sr_reduc    = model_reduction(sr_full, lr_full, m_par)
```
which calculates the unconditional covariance matrix of all state and control variables and rewrites the coefficients of the value functions and the copula as linear combinations of some underlying factors. Only those factors that have eigenvalues above the precision predefined in `sr_full.n_par.compress_critC` (controls, i.e., marginal value functions) and `sr_full.n_par.compress_critS` (states, i.e., the copula) are retained.
!!! warning
    After model reduction, `sr_reduc.indexes_r` contains the indexes that map correctly into the states/controls used in `LOMstate` and `State2Control`.


## Model solution after a parameter change / after reduction
This smaller model (or any model after a parameter change that doesn't affect the steady state) can be solved quickly using a factorization result from [Bayer, Born, and Luetticke (2024, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20201875) running
```
lr_reduc    = update_model(sr_reduc, lr_full, m_par)
```
In the background, this calls [`BASEforHANK.PerturbationSolution.LinearSolution_reduced_system()`](@ref), which only updates the Jacobian entries that regard the **aggregate** model. (Note that both [`BASEforHANK.PerturbationSolution.LinearSolution()`](@ref) and [`BASEforHANK.PerturbationSolution.LinearSolution_reduced_system()`](@ref) call [`BASEforHANK.PerturbationSolution.SolveDiffEq()`](@ref) to obtain a solution to the linearized difference equation.)

This model update step takes about 100ms on a standard computer for the medium size resolution used as a default in the example code.

## Estimation of model parameters
Having obtained `SteadyResults` `sr_reduc` and `LinearResults` `lr_reduc`, the command
```
er_mode = find_mode(sr_reduc, lr_reduc, m_par, e_set)
```
computes the mode of the likelihood, i.e., the parameter vector that maximizes the probability of observing the data given the model, and saves the results in `er_mode`, an instance of `struct` `EstimResults`
(see [`BASEforHANK.Estimation.mode_finding()`](@ref)). We use the Kalman filter to compute the likelihood, and the package `Optim` for optimization. Settings for the estimation can be adjusted in the `struct` [`EstimationSettings`](@ref).

!!! warning
    By default, the flag `estimate_model` in the `struct` [`EstimationSettings`](@ref) is set to `false`. Depending on the computing power available, finding the mode of the likelihood can take several hours to run through. The mode finder might also seem frozen after finishing the optimization but the computation of the Hessian for the large model is involved and can take a long time for the large model. For instructional purposes, we therefore set `e_set.compute_hessian = false` by default and load the Hessian from a save file. For a proper estimation, this has to be set to true. We also save an intermediate step before computing the Hessian in case you are only interested in the mode itself.

Lastly,
```
sample_posterior(sr_reduc, lr_reduc, er_mode, m_par, e_set)
```
uses a Markov Chain Monte Carlo method to trace out the posterior probabilites of the estimated parameters.
The final estimates (and further results) are saved in a file with the name given by the field `save_posterior_file`
in the `struct` `EstimationSettings` (instantiated in `e_set`).

!!! note
    The module `BASEforHANK` creates the estimation settings `e_set` in its main script (when it is initialized),
    so changes to the `struct` `EstimationSettings` are only effective *before* `using BASEforHANK`. Make sure
    that all file paths specified in [`EstimationSettings`](@ref) are correct relative to your script's position.
