"""
    find_steadystate(m_par)

Find the stationary equilibrium capital stock as well as the associated marginal value
functions and the stationary distribution of idiosyncratic states.

This function solves for the market clearing capital stock in the Aiyagari model with
idiosyncratic income risk. That is, it uses [`CustomBrent()`](@ref) to find the root of the
excess capital demand function, which is defined in [`Kdiff()`](@ref). It does so first on a
coarse grid and then on the actual grid.

# Arguments

  - `m_par::ModelParameters`

# Returns

  - `KSS`: Steady-state capital stock
  - `WbSS`, `WkSS`: Marginal value functions
  - `distrSS::Array{Float64,3}`: Steady-state distribution of idiosyncratic states, computed
    by [`Ksupply()`](@ref)
  - `n_par::NumericalParameters`, `m_par::ModelParameters`
"""
function find_steadystate(m_par)

    ## ------------------------------------------------------------------------------------
    ## Step 0: Take care of complete markets case
    ## ------------------------------------------------------------------------------------

    n_par = NumericalParameters(; m_par = m_par)

    if typeof(n_par.model) == CompleteMarkets
        @assert @isdefined(CompMarketsCapital) "Complete Markets Model requires CompMarketsCapital function."

        rSS = (1.0 .- m_par.β) ./ m_par.β  # complete markets interest rate
        KSS = CompMarketsCapital(rSS, m_par)
        WbSS = ones(1, 1, 1)
        WkSS = ones(1, 1, 1)
        distrSS = ones(1, 1, 1)
        n_par = NumericalParameters(;
            m_par = m_par,
            naggrstates = length(state_names),
            naggrcontrols = length(control_names),
            aggr_names = aggr_names,
            distr_names = distr_names,
        )
        return KSS, WbSS, WkSS, distrSS, n_par, m_par
    end

    ## ------------------------------------------------------------------------------------
    ## Step 1: Find the stationary equilibrium for coarse grid
    ## This step is just to get a good starting guess for the actual grid.
    ## ------------------------------------------------------------------------------------

    ## Income process and income grids ----------------------------------------------------

    # Read out numerical parameters for starting guess solution with reduced income grid.
    n_par = NumericalParameters(;
        m_par = m_par,
        ϵ = 1e-6,
        nh = NumericalParameters().nh_coarse,
        nk = NumericalParameters().nk_coarse,
        nb = NumericalParameters().nb_coarse,
    )

    ## Capital stock guesses --------------------------------------------------------------

    if @isdefined(CompMarketsCapital)
        if n_par.verbose
            @printf "CompMarketsCapital function is defined, used for guesses.\n"
        end
        rKmin = n_par.rKmin_coarse
        rKmax = (1.0 .- m_par.β) ./ m_par.β - 0.0025 # complete markets interest rate
        Kmin = CompMarketsCapital(rKmax, m_par)
        Kmax = CompMarketsCapital(rKmin, m_par)
    else
        if n_par.verbose
            @printf "CompMarketsCapital function is not defined, using bounds from n_par.\n"
        end
        Kmin = n_par.Kmin_coarse
        Kmax = n_par.Kmax_coarse
    end

    @assert 0.0 < Kmin < Kmax < Inf "Invalid capital stock bounds."
    @assert !isnan(Kmin) && !isnan(Kmax) "Invalid capital stock bounds."

    if n_par.verbose
        @printf "Kmin: %f, Kmax: %f\n" Kmin Kmax
    end

    ## Excess demand function -------------------------------------------------------------

    # Define the excess demand function based on Kdiff from fcn_kdiff.jl, the keyword
    # arguments allow for a faster solution in CustomBrent because the results of the
    # previous iteration can be used as starting values.
    d(
        K,
        initialize::Bool = true,
        Wb_guess = zeros(1, 1, 1),
        Wk_guess = zeros(1, 1, 1),
        distr_guess = n_par.dist_guess,
    ) = Kdiff(
        K,
        n_par,
        m_par;
        initialize = initialize,
        Wb_guess = Wb_guess,
        Wk_guess = Wk_guess,
        distr_guess = distr_guess,
    )

    ## Find equilibrium capital stock on coarse grid --------------------------------------

    if n_par.verbose
        @printf "Find capital stock, coarse income grid.\n"
    end

    BrentOut = CustomBrent(d, Kmin, Kmax)
    KSS = BrentOut[1]

    if n_par.verbose
        @printf "Capital stock, coarse income grid, is: %f\n" KSS
    end

    ## ------------------------------------------------------------------------------------
    ## Step 2: Find the stationary equilibrium for actual grid
    ## ------------------------------------------------------------------------------------

    ## Update numerical parameters --------------------------------------------------------

    # Read out numerical parameters, update to model equations.
    n_par = NumericalParameters(;
        m_par = m_par,
        naggrstates = length(state_names),
        naggrcontrols = length(control_names),
        aggr_names = aggr_names,
        distr_names = distr_names,
    )

    ## Find equilibrium capital stock on actual grid --------------------------------------

    if n_par.verbose
        @printf "Find capital stock, actual income grid.\n"
    end

    BrentOut = CustomBrent(
        d,
        KSS * (1 - n_par.search_range),
        KSS * (1 + n_par.search_range);
        tol = n_par.ϵ,
    )
    KSS = BrentOut[1]
    WbSS = BrentOut[3][2]
    WkSS = BrentOut[3][3]
    distrSS = BrentOut[3][4]

    if n_par.verbose
        @printf "Capital stock, actual income grid, is: %f\n" KSS
    end

    ## ------------------------------------------------------------------------------------
    ## Return results
    ## ------------------------------------------------------------------------------------

    return KSS, WbSS, WkSS, distrSS, n_par, m_par
end
