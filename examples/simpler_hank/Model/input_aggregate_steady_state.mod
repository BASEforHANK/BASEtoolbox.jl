# Steady state values that are given already:
# - outputs of `Ksupply`, in particular: KSS, BSS
# - all variables in args_hh_prob_names
# - all variables in distr_names

## Shocks
ZSS = m_par.Z

RshockSS = 1.0

## Further assumptions (partly also used in args_hh_prob)
mcSS = 1.0 ./ m_par.μ
mcwSS = 1.0 ./ m_par.μw
πSS = 1.0
πwSS = 1.0

## Variables (that are not already defined)

# nominal interest rates
RBSS = m_par.RRB .* πSS
RLSS = RBSS
RDSS = RRDSS .* πSS

# production side
wFSS = wage(mcSS, ZSS, KSS, NSS, m_par)
YSS = output(ZSS, KSS, NSS, m_par)
ISS = m_par.δ_0 * KSS
Π_FSS = (1.0 - mcSS) .* YSS

# financial market
BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)
BgovSS = BSS
TotalAssetsSS = BSS + qSS * KSS

# fiscal side
RK_before_taxesSS = ((RKSS - 1.0) ./ (1.0 - (TkSS - 1.0))) + 1.0

# jointly determine C, T, G (interacted through consumption tax)
# resource constaint, plugged in government budget constraint and tax revenues
CSS =
    (
        YSS - ISS - (RRDSS .- RRLSS) * BDSS + (RRLSS - 1.0) * BgovSS - (
            (TbarSS .- 1.0) .* (wHSS .* NSS + Π_ESS + Π_USS) +
            (TkSS .- 1.0) .* (RK_before_taxesSS .- 1.0) .* KSS
        )
    ) ./ (1.0 + (TcSS .- 1.0))

# tax revenues
TSS =
    (TbarSS .- 1.0) .* (wHSS .* NSS + Π_ESS + Π_USS) +
    (TcSS .- 1.0) .* CSS +
    (TkSS .- 1.0) .* (RK_before_taxesSS .- 1.0) .* KSS

# government spending from budget constraint
GSS = TSS - (RRLSS - 1.0) * BgovSS

## Lags
BgovlagSS = BgovSS
wFlagSS = wFSS
qlagSS = qSS
TbarlagSS = TbarSS
