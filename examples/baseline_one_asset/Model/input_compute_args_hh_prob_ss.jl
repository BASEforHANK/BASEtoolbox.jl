"""
    compute_args_hh_prob_ss(K, m_par, n_par)

This function calculates the arguments (parameters, variables) that are needed for the
household problem, stored in `args_hh_prob` based on the list `args_hh_prob_names` from
`input_aggregate_names.jl`. These arguments will be packed and unpacked using the macros
[`@write_args_hh_prob`](@ref) (or [`@write_args_hh_prob_ss`](@ref)) and
[`@read_args_hh_prob`](@ref) (or [`@read_args_hh_prob_ss`](@ref)).

This function takes in the current capital stock `K` and the model and numerical parameters
as well as (user-specified) functions and then computes the relevant steady state values.
That means, all the necessary arguments must follow from `K` or steady state assumptions in
the parameter structs.

The function is used in q [`BASEforHANK.SteadyState.Kdiff()`](@ref) to find the stationary
equilibrium of the household block of the model and in
[`BASEforHANK.PerturbationSolution.prepare_linearization()`](@ref) to prepare the
linearization of the model.

# Arguments

  - `K::Float64`: Capital stock
  - `m_par::ModelParameters`, `n_par::NumericalParameters`

# Returns

  - `args_hh_prob::Vector`: Vector of arguments for the household problem, see list
    `args_hh_prob_names` in `input_aggregate_names.jl` and the macros `@write_args_hh_prob`
    and `@write_args_hh_prob_ss`
"""
function compute_args_hh_prob_ss(K, m_par, n_par)

    # Stationary distribution of productivity
    distr_h = (n_par.Π ^ 1000)[1, :]

    # Assumption on markup and markdown
    mc = 1.0 ./ m_par.μ
    mcw = 1.0 ./ m_par.μw

    # Assumption on TFP
    Z = m_par.Z

    # Assumption on Htilde in steady state
    Htilde = n_par.Htilde

    # Assumption on q in steady state
    q = m_par.q

    # Assumption on taxes
    Tlev = m_par.Tlev
    Tprog = m_par.Tprog
    Tc = m_par.Tc

    # Assumption
    σ = m_par.σ

    # Assumption on Hprog in steady state
    Hprog = dot(
        distr_h[1:(end .- 1)],
        (n_par.grid_h[1:(end .- 1)] ./ Htilde) .^ scale_Hprog((Tprog .- 1.0), m_par),
    )
    @assert Hprog ≈ 1.0

    # Calculate market clearing labor supply
    N = find_zero(
        N -> labor_market_clearing_ss(
            N,
            mc,
            mcw,
            Z,
            K,
            Hprog,
            (Tlev .- 1.0),
            (Tprog .- 1.0),
            (Tc .- 1.0),
            m_par,
            m_par.scale_prog,
        ),
        1.0,
    )

    # Calculate interest rate using interest function
    RK = 1.0 .+ interest(mc, Z, K, N, m_par)

    # Calculate wage that firms face using wage function
    wF = wage(mc, Z, K, N, m_par)

    # Calculate output using output function
    Y = output(Z, K, N, m_par)

    # No-arbitrage condition for returns
    RRL = RK

    # Assumption on RRD in steady state
    RRD = borrowing_rate_ss(RRL, m_par)

    # Calculate profits using profits function
    Π_E = profits_E_ss(mc, Y, RRL, m_par)

    # Calculate wage that households face
    wH = mcw .* wF

    # Calculate union profits using union profits function
    Π_U = profits_U(wF, wH, N)

    # Calculate average tax rate using av_labor_tax_rate function
    Tbar =
        1.0 + av_labor_tax_rate(
            n_par,
            m_par,
            wH .* N ./ Hprog,
            (Tlev .- 1.0),
            (Tprog .- 1.0),
            Π_E,
            Htilde,
            distr_h,
        )

    # Package variables for the household block, see list args_hh_prob_names
    @write_args_hh_prob()

    # Return
    return args_hh_prob
end
