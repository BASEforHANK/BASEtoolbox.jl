@doc raw"""
    Kdiff(
        KD,
        n_par,
        m_par,
        initialize = true,
        Wb_guess = zeros(1, 1, 1),
        Wk_guess = zeros(1, 1, 1),
        distr_guess = zeros(1, 1, 1)
    )

This function is used to find the stationary equilibrium of the household block of the
model, in particular, it is used in [`find_steadystate()`](@ref).

Calculate the difference between the capital stock that is demanded/assumed (KD) and the
capital stock that prevails under that demanded capital stock's implied prices when
households face idiosyncratic income risk (Aiyagari model).

Requires global functions from the IncomesETC module and [`Ksupply()`](@ref).

# Arguments
- `KD::Float64`: Assumed capital demand (guess)
- `n_par::NumericalParameters`, `m_par::ModelParameters`
- `initialize::Bool = true`: If true, initialize the marginal value functions and
    stationary distribution, otherwise use the provided guesses that follow. Providing the
    guesses can be used in [`CustomBrent()`](@ref) to speed up the solution since the
    results from the previous iteration can be used as starting values.
- `Wb_guess::AbstractArray = zeros(1, 1, 1)`: Guess for marginal value of liquid assets
- `Wk_guess::AbstractArray = zeros(1, 1, 1)`: Guess for marginal value of illiquid
    assets
- `distr_guess::AbstractArray = zeros(1, 1, 1)`: Guess for stationary distribution

# Returns
- `diff::Float64`: Difference between the demanded and supplied capital stock
- `Wb::AbstractArray`: Marginal value of liquid assets, implied by capital demand
- `Wk::AbstractArray`: Marginal value of illiquid assets, implied by capital demand
- `distr::AbstractArray`: Stationary distribution of idiosyncratic states, implied by
    capital demand
"""
function Kdiff(
    KD::Float64,
    n_par,
    m_par;
    initialize::Bool = true,
    Wb_guess::AbstractArray = zeros(1, 1, 1),
    Wk_guess::AbstractArray = zeros(1, 1, 1),
    distr_guess::AbstractArray = zeros(1, 1, 1),
)

    ## ------------------------------------------------------------------------------------
    ## Step 1: Calculate aggregates implied by capital demand
    ## ------------------------------------------------------------------------------------

    # See module IncomesETC for details on the function compute_args_hh_prob_ss
    args_hh_prob = compute_args_hh_prob_ss(KD, m_par, n_par)

    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Step 2: Calculate incomes of households given aggregates from above
    ## ------------------------------------------------------------------------------------

    # Net incomes of households
    net_income, _, eff_int = incomes(n_par, m_par, args_hh_prob)

    ## ------------------------------------------------------------------------------------
    ## Step 3: Initialize marginal value functions and stationary distribution
    ## ------------------------------------------------------------------------------------

    # If initialize is true, set up guesses here, otherwise use keyword arguments.
    if initialize

        # Guess consumption: consume factor income (if positive)
        net_labor_union_inc_GHH = net_income[1]
        rental_inc = net_income[2]
        liquid_asset_inc = net_income[3]

        # assume that the consumption guess has the following form
        x_guess =
            net_labor_union_inc_GHH .+ rental_inc .* (n_par.mesh_k .* (RK .- 1.0) .> 0) .+
            liquid_asset_inc .* (n_par.mesh_b .> 0)
        if any(any(x_guess .< 0.0))
            @warn "Negative consumption guess. Potentially reduce Kmax_coarse!"
        end

        # Based on consumption guess update the marginal value of liquid assets
        Wb = eff_int .* mutil(x_guess, m_par)

        # Based on consumption guess update the marginal value of illiquid assets
        Wk = (RK + m_par.λ) .* mutil(x_guess, m_par)

        # Guess stationary distribution
        distr = n_par.dist_guess
    else

        # Use provided keyword arguments
        Wb = Wb_guess
        Wk = Wk_guess
        distr = distr_guess
    end

    ## ------------------------------------------------------------------------------------
    ## Step 4: Calculate supply of capital given aggregates from above
    ## ------------------------------------------------------------------------------------

    # Call function Ksupply from fcn_ksupply.jl
    KsupplyOut = Ksupply(args_hh_prob, n_par, m_par, Wb, Wk, distr, net_income, eff_int)

    # Unpack results
    KS = KsupplyOut[1]
    Wb = KsupplyOut[end - 2]
    Wk = KsupplyOut[end - 1]
    distr = KsupplyOut[end]

    # Calculate excess supply of funds
    diff = KS - KD

    ## ------------------------------------------------------------------------------------
    ## Return results
    ## ------------------------------------------------------------------------------------

    return diff, Wb, Wk, distr
end
