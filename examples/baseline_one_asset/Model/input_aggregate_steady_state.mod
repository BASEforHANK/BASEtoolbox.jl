# Steady state values that are given already:
# - outputs of `Ksupply`, in particular: KSS, BSS
# - all variables in args_hh_prob_names
# - all variables in distr_names

## Shocks
ZSS = m_par.Z
ZISS = 1.0
μSS = m_par.μ
μwSS = m_par.μw
ASS = 1.0
RshockSS = 1.0
GshockSS = 1.0
TprogshockSS = 1.0
SshockSS = 1.0

## Growth rates
YgrowthSS = 1.0
BgovgrowthSS = 1.0
IgrowthSS = 1.0
wgrowthSS = 1.0
CgrowthSS = 1.0
TgrowthSS = 1.0

## Further assumptions (partly also used in args_hh_prob)
mcSS = 1.0 ./ μSS
mcwSS = 1.0 ./ μwSS
πSS = 1.0
πwSS = 1.0
uSS = 1.0

## Variables (that are not already defined)

# nominal interest rates
RBSS = RRLSS .* πSS
RLSS = RBSS
RDSS = RRDSS .* πSS

# production side
wFSS = wage(mcSS, ZSS, KSS, NSS, m_par)
YSS = output(ZSS, KSS, NSS, m_par)
ISS = m_par.δ_0 * KSS
Π_FSS = (1.0 - mcSS) .* YSS

# financial market
LPSS = RKSS / (RBSS / πSS)
LPXASS = LPSS
BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)
qΠSS = (m_par.ωΠ .* Π_FSS) ./ (RBSS / πSS .- 1 .+ m_par.ιΠ) + 1.0
BgovSS = BSS .- qΠSS .+ 1.0
TotalAssetsSS = BSS + qSS * KSS

# fiscal side
RK_before_taxesSS = ((RKSS - 1.0) ./ (1.0 - (TkSS - 1.0))) + 1.0
TRSS = transfer_scheme(n_par, m_par, args_hh_prob; distr_h = distr_hSS)

# jointly determine C, T, G (interacted through consumption tax)
# resource constaint, plugged in government budget constraint and tax revenues
CSS =
    (
        YSS - ISS - (RRDSS .- RRLSS) * BDSS + (RRLSS - 1.0) * BgovSS - (
            (TbarSS .- 1.0) .* (wHSS .* NSS + Π_ESS + Π_USS) +
            (TkSS .- 1.0) .* (RK_before_taxesSS .- 1.0) .* KSS - (TRSS .- 1.0)
        )
    ) ./ (1.0 + (TcSS .- 1.0))

# tax revenues
TSS =
    (TbarSS .- 1.0) .* (wHSS .* NSS + Π_ESS + Π_USS) +
    (TcSS .- 1.0) .* CSS +
    (TkSS .- 1.0) .* (RK_before_taxesSS .- 1.0) .* KSS - (TRSS .- 1.0)

# government spending from budget constraint
GSS = TSS - (RRLSS - 1.0) * BgovSS

# remaining aggregates
BYSS = BSS / YSS
TYSS = TSS / YSS

## Lags
YlagSS = YSS
BgovlagSS = BgovSS
TlagSS = TSS
IlagSS = ISS
wFlagSS = wFSS
qlagSS = qSS
ClagSS = CSS
TbarlagSS = TbarSS
TproglagSS = TprogSS
qΠlagSS = qΠSS
