function av_labor_tax_rate(
    n_par,
    m_par,
    labor_compensation,
    τlev,
    τprog,
    Π_E,
    Htilde,
    distr_h,
)

    # Compute aggregate tax base of labor income tax
    tax_base = labor_compensation .+ Π_E

    # Gross labor income of workers, see eq. (Gross income) and (Tax func)
    g_labor_inc = labor_compensation .* (n_par.grid_h / Htilde) .^ scale_Hprog(τprog, m_par)

    # Gross entrepreneurial profits, see eq. (Gross income) and (Tax func)
    g_e_profits = n_par.grid_h[end] .* Π_E

    # Gross income
    g_inc = g_labor_inc
    g_inc[end] = g_e_profits

    # Net income, see eq. (Gross income) and (Tax func)
    n_inc = labor_tax_f.(g_inc, τlev, τprog, tax_base, m_par.scale_prog)

    # Tax revenues
    taxrev = g_inc - n_inc

    # Average tax rate
    τbar = dot(distr_h, taxrev) ./ dot(distr_h, g_inc)

    return τbar
end
