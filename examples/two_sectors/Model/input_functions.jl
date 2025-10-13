## Production -----------------------------------------------------------------------------

# Production function
output(H, S, m_par) = H^m_par.α * S^(1 - m_par.α)

# Real wages that firms pay
wage(mc, Z, H, S, m_par) = mc * (1 - m_par.α) * Z * (H / S)^m_par.α

# Real rental rate of capital, absent utilization-adjusted depreciation
interest(Z, H, S, m_par) = m_par.α * Z * (H / S)^(m_par.α - 1.0) - m_par.δ_0

## Profits --------------------------------------------------------------------------------

# Steady state payout of profits to 'entrepreneur' (housing firm profits are zero in SS)
profits_E_ss(mc, Y, m_par) = (1.0 .- mc) .* Y .* (1 - m_par.α)

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
    S = N .* Z
    H = K .* Z
    wF = wage(mc, Z, H, S, m_par)
    wH = mcw .* wF

    Y = output(H, S, m_par)
    Π_E = profits_E_ss(mc, Y, m_par)
    tax_base = wH .* N ./ Hprog .+ Π_E

    return (N .- labor_supply(wH, Hprog, τlev, τprog, τc, m_par, tax_base, scaling))
end

## Transfers ------------------------------------------------------------------------------

function transfer_scheme(n_par, m_par, args_hh_prob)
    return zeros(size(n_par.grid_h))
end
