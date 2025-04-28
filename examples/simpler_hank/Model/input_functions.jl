## Production -----------------------------------------------------------------------------

# Production function
output(Z, K, N, m_par) = Z * K^m_par.α * N^(1 - m_par.α)

# Real wages that firms pay
wage(mc, Z, K, N, m_par) = mc * (1 - m_par.α) * Z * (K / N)^m_par.α

# Real rental rate of capital, absent utilization-adjusted depreciation
interest(mc, Z, K, N, m_par) = mc * m_par.α * Z * (K / N)^(m_par.α - 1.0) - m_par.δ_0

## Profits --------------------------------------------------------------------------------

# Steady state payout to entrepreneurs
profits_E_ss(mc, Y) = (1.0 .- mc) .* Y

## Union ----------------------------------------------------------------------------------

# Union profits
profits_U(wF, wH, N) = (wF .- wH) .* N

## Financial markets ----------------------------------------------------------------------

# Borrowing rate, as function of lending rate and parameters
borrowing_rate_ss(RRL, m_par) = RRL .+ m_par.Rbar

## Labor market clearing ------------------------------------------------------------------

# Labor market clearing, solving for equilibrium labor given the firm side: this could also
# derived closed form and adjusted in `compute_args_hh_prob_ss()` accordingly.
function labor_market_clearing_ss(
    N,
    mc,
    mcw,
    Z,
    K,
    Hprog,
    τlev,
    τprog,
    τc,
    m_par,
    scaling::Bool,
)
    wF = wage(mc, Z, K, N, m_par)
    wH = mcw .* wF

    Y = output(Z, K, N, m_par)
    Π_E = profits_E_ss(mc, Y)
    tax_base = wH .* N ./ Hprog .+ Π_E

    return (N .- labor_supply(wH, Hprog, τlev, τprog, τc, m_par, tax_base, scaling))
end

## Optional functions ---------------------------------------------------------------------

# CompMarketsCapital is used in find_steadystate.jl to improve initial guesses, if not
# provided, the model will use the initial guess for K from the model parameters. First, get
# K/N (capital intensity) from rearranging interest function, then, compute full insurance
# labor supply and return K given these assumptions. This is hard-coded and should be
# adjusted if the model changes.
function CompMarketsCapital(rK, m_par)
    Z = m_par.Z
    mc = 1 / m_par.μ
    mcw = 1 / m_par.μw
    K_over_N = ((rK + m_par.δ_0) / (m_par.α * Z * mc))^(1 / (m_par.α - 1))
    wF = wage(mc, Z, K_over_N, 1.0, m_par)
    wH = mcw .* wF
    Hprog = 1.0
    N = labor_supply(
        wH,
        Hprog,
        (m_par.Tlev .- 1.0),
        (m_par.Tprog .- 1.0),
        (m_par.Tc .- 1.0),
        m_par,
        0.0,
        false,
    )
    return K_over_N * N
end
