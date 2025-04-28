#=

Template function for prepare_linearization.jl

Explanation: During the preprocessing step, `PreprocessInputs.jl` reads in this template
file and copies the content from `input_aggregate_steady_state.mod` into the code below at
the line marked with "# aggregate steady state marker". The code block is then written to
`prepare_linearization_generated.jl` in the `generated_fcns` directory.

=#

"""
    prepare_linearization(K, Wb, Wk, distr, n_par, m_par)

Given the stationary equilibrium of the household side, computed in
[`find_steadystate()`](@ref), this function performs several steps:

  - Step 1: compute the stationary equilibrium.
  - Step 2: perform the dimensionality reduction of the marginal value functions as well as
    the distribution.
  - Step 3: compute the aggregate steady state from `input_aggregate_steady_state.mod`.
  - Step 4: produce indexes to access the variables in the linearized model.
  - Step 5: return the results.

# Arguments

  - `K::Float64`: steady-state capital stock
  - `Wb::Array{Float64,3}`, `Wk::Array{Float64,3}`: steady-state marginal value functions
  - `distr::Array{Float64,3}`: steady-state distribution of idiosyncratic states
  - `n_par::NumericalParameters`,`m_par::ModelParameters`

# Returns

  - `XSS::Array{Float64,1}`, `XSSaggr::Array{Float64,1}`: steady state vectors produced by
    [`@writeXSS()`](@ref)
  - `indexes`, `indexes_aggr`: `struct`s for accessing `XSS`,`XSSaggr` by variable names,
    produced by [`@make_fn()`](@ref), [`@make_fnaggr()`](@ref)
  - `compressionIndexes::Array{Array{Int,1},1}`: indexes for compressed marginal value
    functions (`V_m` and `V_k`)
  - `n_par::NumericalParameters`, `m_par::ModelParameters`: updated parameters
  - `CDFSS`, `CDF_bSS`, `CDF_kSS`, `CDF_hSS`: cumulative distribution functions (joint and
    marginals)
  - `distrSS::Array{Float64,3}`: steady state distribution of idiosyncratic states, computed
    by [`Ksupply()`](@ref)
"""
function prepare_linearization(
    K::Float64,
    Wb::Array{Float64,3},
    Wk::Array{Float64,3},
    distr::Array{Float64,3},
    n_par::NumericalParameters,
    m_par::ModelParameters,
)

    # Guarantee stability of the scope of the function
    KSS = copy(K)
    WbSS = copy(Wb)
    WkSS = copy(Wk)
    distrSS = copy(distr)

    ## ------------------------------------------------------------------------------------
    ## Step 1: Evaluate stationary equilibrium to calculate steady state variable values
    ## ------------------------------------------------------------------------------------

    ## Aggregate part ---------------------------------------------------------------------

    # See module IncomesETC for details on the function compute_args_hh_prob_ss
    args_hh_prob = compute_args_hh_prob_ss(KSS, m_par, n_par)

    @read_args_hh_prob_ss()

    # Store the arguments in a named tuple, for checking consistency
    args_hh_prob_tuple = NamedTuple{Tuple(Symbol.(args_hh_prob_names))}(args_hh_prob)

    ## Incomes ----------------------------------------------------------------------------

    # Net incomes of households
    net_income, gross_income, eff_int = incomes(n_par, m_par, args_hh_prob)

    ## Idiosyncratic part -----------------------------------------------------------------

    # Solution to household problem in the stationary equilibrium, given args_hh_prob
    KSS,
    BSS,
    ΓSS,
    Γ_aSS,
    Γ_nSS,
    x_a_starSS,
    b_a_starSS,
    k_a_starSS,
    x_n_starSS,
    b_n_starSS,
    WbSS,
    WkSS,
    distrSS = Ksupply(args_hh_prob, n_par, m_par, WbSS, WkSS, distrSS, net_income, eff_int)

    WbSS .*= eff_int

    # Distributional summary statistics
    distr_bSS,
    distr_kSS,
    distr_hSS,
    TOP10WshareSS,
    TOP10IshareSS,
    TOP10InetshareSS,
    GiniWSS,
    GiniCSS,
    sdlogySS = distrSummaries(
        distrSS,
        qSS,
        x_a_starSS,
        x_n_starSS,
        n_par,
        net_income,
        gross_income,
        m_par,
    )

    ## ------------------------------------------------------------------------------------
    ## Step 2: Dimensionality reduction
    ## ------------------------------------------------------------------------------------

    ## DCT coefficients for marginal value functions --------------------------------------

    # Identify the co-linear variables in arguments of household problem
    exclude_list = ["N", "Hprog", "Htilde", "RRD"]
    include_list_idx = findall(x -> x ∉ exclude_list, args_hh_prob_names)

    if typeof(n_par.model) == CompleteMarkets
        indb = [1 2]
        indk = [1 2]
    else
        # Perform first stage reduction
        indk, indb, _ = first_stage_reduction(
            WkSS,
            WbSS,
            Γ_aSS,
            Γ_nSS,
            b_a_starSS,
            k_a_starSS,
            b_n_starSS,
            args_hh_prob,
            include_list_idx,
            n_par,
            m_par,
        )
    end

    WbSS = log.(invmutil(WbSS, m_par))
    WkSS = log.(invmutil(WkSS, m_par))
    compressionIndexesWb =
        typeof(n_par.model) == CompleteMarkets ? [] : sort(unique(vcat(indb...)))
    compressionIndexesWk =
        typeof(n_par.model) == TwoAsset ? sort(unique(vcat(indk...))) : []

    ## Polynomials for copula perturbation ------------------------------------------------

    SELECT = [
        (!((i == 1) & (j == 1)) & !((k == 1) & (j == 1)) & !((k == 1) & (i == 1))) for
        i = 1:(n_par.nb_copula), j = 1:(n_par.nk_copula), k = 1:(n_par.nh_copula)
    ]

    # store indices of selected coeffs
    compressionIndexesCOP = findall(SELECT[:])

    ## Store compression indexes ----------------------------------------------------------

    # Container to store all retained coefficients in one array
    compressionIndexes = Array{Array{Int,1},1}(undef, 3)
    compressionIndexes[1] = compressionIndexesWb
    compressionIndexes[2] = compressionIndexesWk
    compressionIndexes[3] = compressionIndexesCOP

    ## Produce marginals ------------------------------------------------------------------

    # Calculate CDF from PDF
    CDFSS = cumsum(cumsum(cumsum(distrSS; dims = 1); dims = 2); dims = 3)

    # Marginal distribution (pdf) of liquid assets
    distr_bSS = sum(distrSS; dims = (2, 3))[:]

    # Marginal distribution (pdf) of illiquid assets
    distr_kSS = sum(distrSS; dims = (1, 3))[:]

    # Marginal distribution (pdf) of productivity
    distr_hSS = sum(distrSS; dims = (1, 2))[:]

    # Marginal distribution (cdf) of liquid assets
    CDF_bSS = cumsum(distr_bSS[:])

    # Marginal distribution (cdf) of illiquid assets
    CDF_kSS = cumsum(distr_kSS[:])

    # Marginal distribution (cdf) of productivity
    CDF_hSS = cumsum(distr_hSS[:])

    # Calculate interpolation nodes for the copula as those elements of the marginal
    # distribution that yield close to equal aggregate shares in liquid wealth, illiquid
    # wealth and productivity. Entrepreneur state treated separately.
    @set! n_par.copula_marginal_b =
        n_par.nb == 1 ? distr_bSS :
        copula_marg_equi(distr_bSS, n_par.grid_b, n_par.nb_copula)
    @set! n_par.copula_marginal_k =
        n_par.nk == 1 ? distr_kSS :
        copula_marg_equi(distr_kSS, n_par.grid_k, n_par.nk_copula)
    @set! n_par.copula_marginal_h =
        n_par.nh == 2 ? distr_hSS :
        copula_marg_equi_h(distr_hSS, n_par.grid_h, n_par.nh_copula)

    ## ------------------------------------------------------------------------------------
    ## Step 3: Get the aggregate steady state (`input_aggregate_steady_state.mod`)
    ## ------------------------------------------------------------------------------------

    # DO NOT DELETE OR EDIT NEXT LINE! This is needed for parser.

    # aggregate steady state marker

    # Write steady state values into XSS vector
    @writeXSS

    ## ------------------------------------------------------------------------------------
    ## Step 4: Produce indexes to access the variables in the linearized model
    ## ------------------------------------------------------------------------------------

    # produce indexes to access XSS etc.
    indexes = produce_indexes(
        n_par,
        compressionIndexesWb,
        compressionIndexesWk,
        compressionIndexesCOP,
    )
    indexes_aggr = produce_indexes_aggr(n_par)

    @set! n_par.ntotal =
        length(vcat(compressionIndexes...)) +
        (n_par.nh + n_par.nb + n_par.nk - 3 + n_par.naggr)
    @set! n_par.nstates =
        n_par.nh + n_par.nk + n_par.nb - 3 +
        n_par.naggrstates +
        length(compressionIndexes[3])
    @set! n_par.ncontrols = length(vcat(compressionIndexes[1:2]...)) + n_par.naggrcontrols
    @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
    @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
    @set! n_par.nstates_r = copy(n_par.nstates)
    @set! n_par.ncontrols_r = copy(n_par.ncontrols)
    @set! n_par.ntotal_r = copy(n_par.ntotal)
    @set! n_par.PRightStates = Diagonal(ones(n_par.nstates))
    @set! n_par.PRightAll = Diagonal(ones(n_par.ntotal))

    if n_par.n_agg_eqn != n_par.naggr - length(n_par.distr_names)
        @warn "Inconsistency in number of aggregate variables and equations!"
    end

    ## ------------------------------------------------------------------------------------
    ## Step 5: Check consistency
    ## ------------------------------------------------------------------------------------

    # Check that the steady state values in the aggregate model are consistent
    for (key, value) in pairs(args_hh_prob_tuple)
        if !isapprox(exp(XSS[getfield(indexes, Symbol(key, "SS"))]), value)
            @warn "Inconsistency detected for $key: expected $value, got $(exp(XSS[getfield(indexes, Symbol(key, "SS"))]))"
        end
    end

    ## ------------------------------------------------------------------------------------
    ## Step 6: Return results
    ## ------------------------------------------------------------------------------------

    return XSS,
    XSSaggr,
    indexes,
    indexes_aggr,
    compressionIndexes,
    n_par,
    m_par,
    CDFSS,
    CDF_bSS,
    CDF_kSS,
    CDF_hSS,
    distrSS
end

## ----------------------------------------------------------------------------------------
## Helper functions
## ----------------------------------------------------------------------------------------

"""
    copula_marg_equi_h(distr_i, grid_i, nx)
"""
function copula_marg_equi_h(distr_i, grid_i, nx)
    CDF_i = cumsum(distr_i[:])
    aux_marginal = collect(range(CDF_i[1]; stop = CDF_i[end], length = nx))

    x2 = 1.0
    for i = 2:(nx - 1)
        equi(x1) = equishares(x1, x2, grid_i[1:(end - 1)], distr_i[1:(end - 1)], nx - 1)
        x2 = find_zero(equi, (1e-9, x2))
        aux_marginal[end - i] = x2
    end

    aux_marginal[end] = CDF_i[end]
    aux_marginal[1] = CDF_i[1]
    aux_marginal[end - 1] = CDF_i[end - 1]
    copula_marginal = copy(aux_marginal)
    jlast = nx - 1
    for i = (nx - 2):-1:1
        j = locate(aux_marginal[i], CDF_i) + 1
        if jlast == j
            j -= 1
        end
        jlast = j
        copula_marginal[i] = CDF_i[j]
    end
    return copula_marginal
end

"""
    copula_marg_equi(distr_i, grid_i, nx)
"""
function copula_marg_equi(distr_i, grid_i, nx)
    CDF_i = cumsum(distr_i[:])
    aux_marginal = collect(range(CDF_i[1]; stop = CDF_i[end], length = nx))

    x2 = 1.0
    for i = 1:(nx - 1)
        equi(x1) = equishares(x1, x2, grid_i, distr_i, nx)
        x2 = find_zero(equi, (eps(), x2))
        aux_marginal[end - i] = x2
    end

    aux_marginal[end] = CDF_i[end]
    aux_marginal[1] = CDF_i[1]
    copula_marginal = copy(aux_marginal)
    jlast = nx
    for i = (nx - 1):-1:1
        j = locate(aux_marginal[i], CDF_i) + 1
        if jlast == j
            j -= 1
        end
        jlast = j
        copula_marginal[i] = CDF_i[j]
    end
    return copula_marginal
end

"""
    equishares(x1, x2, grid_i, distr_i, nx)
"""
function equishares(x1, x2, grid_i, distr_i, nx)
    FN_Wshares = cumsum(grid_i .* distr_i) ./ sum(grid_i .* distr_i)
    Wshares = diff(mylinearinterpolate(cumsum(distr_i), FN_Wshares, [x1; x2]))
    dev_equi = Wshares .- 1.0 ./ nx

    return dev_equi
end
