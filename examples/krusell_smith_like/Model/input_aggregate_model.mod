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

# Constant markups (like RBC model)
mcw = 1.0 ./ m_par.μw
mc = 1.0 ./ m_par.μ

## ----------------------------------------------------------------------------------------
## Aggregate equations
## ----------------------------------------------------------------------------------------

## Fiscal policy --------------------------------------------------------------------------

# Average tax rate, see equation 34 in BBL (here simplified)
F[indexes.Tbar] = (log(Tbar)) - (XSS[indexes.TbarSS])
# This variable needs to be set for the package!

# Progressivity of labor tax, see equation 33a in BBL (here simplified)
F[indexes.Tprog] = (log(Tprog)) - (XSS[indexes.TprogSS])
# This variable needs to be set for the package!

# Level of labor tax, see equation 35 in BBL (typos!), this determines Tlev
F[indexes.Tlev] = (log(Tlev)) - (XSS[indexes.TlevSS])
# This variable needs to be set for the package!

# VAT rate (gross)
F[indexes.Tc] = (log(Tc)) - (XSS[indexes.TcSS])
# This variable needs to be set for the package!

## Labor market ---------------------------------------------------------------------------

# Idiosyncratic income risk (contemporaneous reaction to business cycle)
F[indexes.σ] = (log(σ)) - (XSS[indexes.σSS])
# This variable needs to be set for the package!

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
F[indexes.Hprog] = (log(Hprog)) - (XSS[indexes.HprogSS])
# This variable needs to be set for the package!

## Production -----------------------------------------------------------------------------

# Rate of return on capital
F[indexes.RK] = (log(RK)) - (log(1 + interest(mc, Z, K, N, m_par)))
# This variable needs to be set for the package!

# Wages that firms pay
F[indexes.wF] = (log(wF)) - (log(wage(mc, Z, K, N, m_par)))

# Firm profits
F[indexes.Π_F] = (log(Π_F)) - (log(Y * (1.0 - mc) + (KPrime - (1.0 - m_par.δ_0) * K) - I))

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

# Real rates on liquid assets
F[indexes.RRL] = (log(RRL)) - (log(RK))
# This variable needs to be set for the package!

# Real rates on liquid debt
F[indexes.RRD] = (log(RRD)) - (log(borrowing_rate_ss(RRL, m_par)))
# This variable needs to be set for the package!

## Market clearing ------------------------------------------------------------------------

# Resource constraint
F[indexes.C] = (log(Y - I)) - (log(C))

## Additional definitions -----------------------------------------------------------------

# Capital by accounting identity
F[indexes.K] = (log(K)) - (log(TotalAssets)) # qlag * K = TotalAssets (but q = 1 always)

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

# Total assets market clearing
F[indexes.TotalAssets] = (log(TotalAssets)) - (XSS[indexes.TotalAssetsSS])

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
