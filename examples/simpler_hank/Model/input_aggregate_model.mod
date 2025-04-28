#=

This file contains the aggregate model equations. That is, everything but the household
planning problem which is described by one EGM backward step and one forward iteration of
the distribution.

Model equations take the form F[equation number] = (lhs) - (rhs)

Equation numbers are generated automatically and stored in the index struct. For this the
corresponding variable needs to be in the list of states or controls.

=#

## ----------------------------------------------------------------------------------------
## Auxiliary variables
## ----------------------------------------------------------------------------------------

## Taxation -------------------------------------------------------------------------------

# Mass of households in each productivity state, distribution is (nb, nk, nh)
distr_h = sum(distrSS; dims = (1, 2))

## ----------------------------------------------------------------------------------------
## Aggregate equations
## ----------------------------------------------------------------------------------------

## Lagged variables -----------------------------------------------------------------------

F[indexes.Bgovlag] = (log(BgovlagPrime)) - (log(Bgov))
F[indexes.wFlag] = (log(wFlagPrime)) - (log(wF))
F[indexes.qlag] = (log(qlagPrime)) - (log(q))
F[indexes.Tbarlag] = (log(TbarlagPrime)) - (log(Tbar))

## Fiscal policy --------------------------------------------------------------------------

# Deficit rule, see equation 33 in BBL
F[indexes.π] =
    (log(BgovPrime / Bgov)) - (
        -m_par.γ_B * (log(Bgov) - XSS[indexes.BgovSS]) +
        m_par.γ_Y * log(Y / exp(XSS[indexes.YSS])) +
        m_par.γ_π * log(π)
    )

# Average tax rate, see equation 34 in BBL (here simplified)
F[indexes.Tbar] =
    (log(Tbar .- 1.0)) - (
        m_par.ρ_τ * log(Tbarlag .- 1.0) +
        (1.0 - m_par.ρ_τ) * (log(exp(XSS[indexes.TbarSS]) .- 1.0)) +
        (1.0 - m_par.ρ_τ) * m_par.γ_Yτ * log(Y / exp(XSS[indexes.YSS])) +
        (1.0 - m_par.ρ_τ) * m_par.γ_Bτ * (log(Bgov) - log(Bgovlag))
    )
# This variable needs to be set for the package!

# Progressivity of labor tax, see equation 33a in BBL (here simplified)
F[indexes.Tprog] = (log(Tprog)) - (XSS[indexes.TprogSS])
# This variable needs to be set for the package!

# Level of labor tax, see equation 35 in BBL (typos!), this determines Tlev
F[indexes.Tlev] =
    (Tbar .- 1.0) - (av_labor_tax_rate(
        n_par,
        m_par,
        wH * N / Hprog,
        (Tlev .- 1.0),
        (Tprog .- 1.0),
        Π_E,
        Htilde,
        distr_h,
    ))
# This variable needs to be set for the package!

# Government budget constraint, see below equation 35 in BBL
F[indexes.G] = (log(G)) - (log(BgovPrime + T - RB / π * Bgov))

# Total goverment tax revenues, see below equation 35 in BBL
F[indexes.T] =
    (log(T)) - (log((Tbar .- 1.0) * (wH * N) + (Tbar .- 1.0) * Π_E + (Tbar .- 1.0) * Π_U))

# VAT rate (gross)
F[indexes.Tc] = (log(Tc)) - (XSS[indexes.TcSS])
# This variable needs to be set for the package!

## Monetary policy ------------------------------------------------------------------------

# Taylor rule, see equation 32 in BBL (here simplified)
F[indexes.RB] =
    (log(RBPrime)) - (
        XSS[indexes.RBSS] +
        ((1 - m_par.ρ_R) * m_par.θ_π) * log(π) +
        ((1 - m_par.ρ_R) * m_par.θ_Y) * log(Y / exp(XSS[indexes.YSS])) +
        m_par.ρ_R * (log(RB) - XSS[indexes.RBSS]) +
        log(Rshock)
    )

# Monetary policy shock
F[indexes.Rshock] = (log(RshockPrime)) - (m_par.ρ_Rshock * log(Rshock))

## Labor market ---------------------------------------------------------------------------

# Idiosyncratic income risk (contemporaneous reaction to business cycle)
F[indexes.σ] = (log(σ)) - (XSS[indexes.σSS])
# This variable needs to be set for the package!

# Wage Phillips Curve
F[indexes.mcw] =
    (log(πw) - XSS[indexes.πwSS]) - (
        m_par.κw * (mcw - 1.0 / m_par.μw) +
        m_par.β * ((log(πwPrime) - XSS[indexes.πwSS]) * (NPrime * wFPrime) / (N * wF))
    )

# Definition of real wage inflation
F[indexes.πw] = (log(wF / wFlag)) - (log(πw / π))

# Wages that households receive
F[indexes.wH] = (log(wH)) - (log(mcw * wF))
# This variable needs to be set for the package!

# Union profits
F[indexes.Π_U] = (log(Π_U)) - (log(profits_U(wF, wH, N)))
# This variable needs to be set for the package!

# Labor supply
F[indexes.N] =
    (log(N)) - (log(
        labor_supply(
            wH,
            Hprog,
            (Tlev .- 1.0),
            (Tprog .- 1.0),
            (Tc .- 1.0),
            m_par,
            wH * N / Hprog + Π_E,
            m_par.scale_prog,
        ),
    ))
# This variable needs to be set for the package!

# Hours-weighted average labor productivity, normalized, see equation 19b in BBL
F[indexes.Hprog] =
    (log(Hprog)) - (log(
        dot(
            distr_h[1:(end - 1)],
            (n_par.grid_h[1:(end - 1)] / Htilde) .^ scale_Hprog((Tprog .- 1.0), m_par),
        ),
    ))
# This variable needs to be set for the package!

## Production -----------------------------------------------------------------------------

# Price Phillips Curve
F[indexes.mc] =
    (log(π) - XSS[indexes.πSS]) - (
        m_par.κ * (mc - 1 / m_par.μ) +
        m_par.β * ((log(πPrime) - XSS[indexes.πSS]) * YPrime / Y)
    )

# Rate of return on capital
F[indexes.RK] =
    (log(RK)) - (log(1.0 + interest(mc, Z, K, N, m_par) + m_par.δ_0 - q * m_par.δ_0))
# This variable needs to be set for the package!

# Wages that firms pay
F[indexes.wF] = (log(wF)) - (log(wage(mc, Z, K, N, m_par)))

# Firm profits
F[indexes.Π_F] =
    (log(Π_F)) - (log(Y * (1.0 - mc) + q * (KPrime - (1.0 - m_par.δ_0) * K) - I))

# Distributed profits to entrepreneurs
F[indexes.Π_E] = (log(Π_E)) - (log(Π_F))
# This variable needs to be set for the package!

# Price of capital investment
F[indexes.q] = (log(q)) - (log(m_par.q))
# This variable needs to be set for the package!

# Capital accumulation equation
F[indexes.I] = (log(KPrime)) - (log(K * (1.0 - m_par.δ_0) + I))

# Production function
F[indexes.Y] = (log(Y)) - (log(output(Z, K, N, m_par)))

# TFP
F[indexes.Z] = (log(ZPrime)) - (m_par.ρ_Z * log(Z))

## Asset markets --------------------------------------------------------------------------

# Return on liquid assets
F[indexes.RL] = (log(RL)) - (log(RB))

# Return on liquid debt
F[indexes.RD] = (log(RD)) - (log(RRD .* π))

# Total liquidity demand
F[indexes.Bgov] = (log(B)) - (log(Bgov))

# Real rates on liquid assets
F[indexes.RRL] = (log(RRL)) - (log(RL / π))
# This variable needs to be set for the package!

# Real rates on liquid debt
F[indexes.RRD] = (log(RRD)) - (log(borrowing_rate_ss(RRL, m_par)))
# This variable needs to be set for the package!

## Market clearing ------------------------------------------------------------------------

# Resource constraint
F[indexes.C] = (log(Y - G - I - BD * (RRD .- RRL))) - (log(C))

## Additional definitions -----------------------------------------------------------------

# Total assets by accounting identity
F[indexes.TotalAssets] = (log(TotalAssets)) - (log(qlag * K + B))

## ----------------------------------------------------------------------------------------
## Closing the aggregate model, see documentation for details
## ----------------------------------------------------------------------------------------

#=

Do not delete the following lines of code!

These equations are overwritten in FSYS by the corresponding aggregation equations of
households' decisions. Here, they are simply set to close the aggregate model. This is a
trick that is exploited in the estimation when only the derivatives with respect to
aggregates is needed. These derivatives are still correct since the left-hand-side of the
equations are the same in both the purely aggregate as well as the complete model.

=#

# Scaling factor for individual productivity
F[indexes.Htilde] = (log(Htilde)) - (XSS[indexes.HtildeSS])
# This variable needs to be set for the package!

# Capital market clearing
F[indexes.K] = (log(K)) - (XSS[indexes.KSS])

# Bond market clearing
F[indexes.B] = (log(B)) - (XSS[indexes.BSS])

# IOUs
F[indexes.BD] = (log(BD)) - (XSS[indexes.BDSS])

## ----------------------------------------------------------------------------------------
## Other distributional statistics
## ----------------------------------------------------------------------------------------

#=

# TO BE EXPLAINED! CURRENTLY, THESE ARE HARD-CODED IN FSYS, FIX THIS!

=#

# other distributional statistics not used in other aggregate equations and not changing
# with parameters, but potentially with other aggregate variables are NOT included here.
# They are found in FSYS.
