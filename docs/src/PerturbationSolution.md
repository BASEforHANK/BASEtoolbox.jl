# Linear perturbation around steady state

!!! note
    The main functions of this section are in the Submodule `PerturbationSolution`.

The model is linearized with respect to aggregate variables. For this,
we write the equilibrium conditions in the form of
``F(X,X')=0``, where ``X`` and ``X'`` are (expected) deviations from steady state
in two successive periods. Applying the total differential yields
``A*X' = - B*X``, where ``A``,``B`` are the first derivatives of ``F`` with respect
to ``X'``,``X``. In the standard setting, we use the generalized Schur decomposition [^Klein]
to transform this equation into a linearized observation equation ``d = gx*k`` and
a linearized state transition equation ``k' = hx*k``, where ``k`` is a vector of the
*state* variables and ``d`` is a vector of the *control* variables (``X = \begin{bmatrix} k \\ d \end{bmatrix}``).

In our code, ``F`` is implemented as [`BASEforHANK.PerturbationSolution.Fsys()`](@ref), while differentiating and
solving for ``gx`` and ``hx`` is done in [`BASEforHANK.PerturbationSolution.LinearSolution()`](@ref), called by [`linearize_full_model()`](@ref) returns the results as a `struct` `LinearResults`:

## `linearize_full_model()`
```@docs
linearize_full_model
BASEforHANK.LinearSolution
```
The function `linearize_full_model()` calls `LinearSolution` and stores the results in a `LinearResults` `struct`.

### `LinearSolution()` 
The `LinearSolution` function executes the following steps:

- generate devices to retrieve distribution and marginal value functions from
    compressed states/controls (`Γ` and `DC`,`IDC`)
- calculate the first derivative of [`BASEforHANK.PerturbationSolution.Fsys()`](@ref) with respect to `X` and `XPrime`.
    We use automatic differentiation (implemented in Julia by the package `ForwardDiff`).
    Partial derivatives are calculated using the `ForwardDiff.jacobian()` function.
    We exploit that some partial derivatives have known values (contemporaneous marginal value
    functions and the future marginal distributions) and set them directly instead of calculating them [^BL].

- compute linear observation and state transition equations using the [`BASEforHANK.PerturbationSolution.SolveDiffEq()`](@ref) function

### `SolveDiffEq()`
```@docs
BASEforHANK.PerturbationSolution.SolveDiffEq
```
- apply the model reduction by pre- and post-multiplying the reduction matrix ``\mathcal{P}`` to the Jacobians `A` and `B` that are inputs to [`BASEforHANK.PerturbationSolution.SolveDiffEq()`](@ref). ``\mathcal{P}`` is calculated in [`model_reduction()`](@ref) and stored in `n_par.PRightAll`.
- compute linear observation and state transition equations. The solution algorithm is set
    in `n_par.sol_algo`, with the options `:schur` (mentioned above) and `:litx` [^lit]. The results are matrices that map contemporaneous states to controls [`gx`],
    or contemporaneous states to future states [`hx`]


### `Fsys()`
```@docs
BASEforHANK.PerturbationSolution.Fsys
```
The function [`BASEforHANK.PerturbationSolution.Fsys()`](@ref) proceeds in the following way:
1. set up vector `F`, that contains the errors to all equilibrium conditions. There are as many conditions
    as deviations from steady state (length of `X`,`XPrime`), and conditions are indexed with
    respective model variable in `IndexStruct` `indexes`
2. generate locally all aggregate variables (for both periods) using [`BASEforHANK.Parsing.@generate_equations`](@ref)
3. construct the full-grid marginal distributions, marginal value functions, and the copula
    from the steady-state values and the (compressed) deviations (for the copula, the selection of DCT
    coefficients that can be perturbed ensures that also the perturbed function is a copula)
4. write all equilibrium condition-errors with respect to *aggregate* variables to `F`, using
    [`BASEforHANK.PerturbationSolution.Fsys_agg()`](@ref)
5. compute optimal policies with [`BASEforHANK.SteadyState.EGM_policyupdate()`](@ref), given
    future marginal value functions, prices, and individual incomes. Infer present marginal
    value functions from them (envelope theorem) and set the difference to assumed present
    marginal value functions (in terms of their compressed deviation from steady state)
    as equilibrium condition-errors (*backward iteration of the value function*)
6. compute future marginal distributions and the copula (on the copula grid) from previous distribution and optimal asset policies. Interpolate when necessary. Set difference to assumed future marginal distributions and copula values on the copula nodes as equilibrium condition-errors (*forward iteration of the distribution*)
7. compute distribution summary statistics with [`BASEforHANK.Tools.distrSummaries()`](@ref) and write
    equilibrium conditions with their respective (control) variables
8. return `F`

Note that the copula is treated as the sum of two interpolants. An interpolant based on the steady-state distribution using the full steady-state marginals as a grid and a "deviations"-function that is defined on the copula grid generated in `prepare_linearization()`. The actual interpolation is carried out with [`BASEforHANK.Tools.myinterpolate3()`](@ref). Default setting is trilinear interpolation, the code also allows for 3d-Akima interpolation.

## `model_reduction()`
```@docs
model_reduction
```
The function [`model_reduction()`](@ref) derives the approximate factor representation from a first solution of the heterogeneous agent model.[^BBL] It then stores the matrices that allow to map the factors to the full set of state and control variables. For deriving the factor representation, the function calculates the long run variance-covariance matrix of all states of the model (given its first-stage reduction).


## `update_model()`
```@docs
update_model
BASEforHANK.PerturbationSolution.Fsys_agg
```
The function [`update_model()`](@ref) solves the aggregate model without updating the derivatives of the household/idiosyncratic part. For this purpose the derivatives of [`BASEforHANK.PerturbationSolution.Fsys_agg()`](@ref) are calculated instead of `Fsys()`. This substantially speeds up the solution after a parameter change that only affects aggregates.[^BBL] In particular, if the model is reduced to its approximate factor representation (see above), this generates significant speed gains.


[^Klein]:
    See the paper [Using the generalized Schur form to solve a multivariate linear rational expectations model](https://www.sciencedirect.com/science/article/pii/S0165188999000457) by Paul Klein (JEDC 2000)

[^BL]:
    Contemporaneous marginal value functions are irrelevant for optimal decisions, so
    its effect on other model variables is 0. Due to a rich enough set of prices, the future distribution
    directly only affects the Fokker-Planck equation. For details, see the paper
    [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](https://doi.org/10.3982/QE1243), *Quantitative Economics*, Vol.11(4), November 2020, p. 1253-1288.

[^BBL]:
    For details, see the paper [Shocks, Frictions, and Inequality in US Business Cycles](https://www.benjaminborn.de/files/BBL_Inequality_Sep2023.pdf), *American Economic Review*, forthcoming.

[^lit]:
    Invoking the Implicit Function Theorem, there exist functions ``g`` and ``h`` such that
    ``F\left(\begin{pmatrix} k \\ g(k) \end{pmatrix},\begin{pmatrix} h(k) \\ g(h(k)) \end{pmatrix}\right)=0``.
    Totally differentiating by ``k`` yields ``B \begin{pmatrix}\mathbb{I}\\ Dg \end{pmatrix}+A \begin{pmatrix}\mathbb{I}\\ Dg \end{pmatrix} Dh = 0``. The `:lit`-algorithm solves this equation for ``Dg`` and ``Dh`` iteratively.

    