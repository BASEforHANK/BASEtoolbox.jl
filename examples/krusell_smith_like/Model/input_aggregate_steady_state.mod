# Steady state values that are given already:
# - outputs of `Ksupply`, in particular: KSS, BSS
# - all variables in args_hh_prob_names
# - all variables in distr_names

## Shocks
ZSS = m_par.Z

## Further assumptions (partly also used in args_hh_prob)
mcSS = 1.0 ./ m_par.μ
mcwSS = 1.0 ./ m_par.μw

# production side
wFSS = wage(mcSS, ZSS, KSS, NSS, m_par)
YSS = output(ZSS, KSS, NSS, m_par)
ISS = m_par.δ_0 * KSS
Π_FSS = (1.0 - mcSS) .* YSS

# financial market
BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)
TotalAssetsSS = KSS

# fiscal side
RK_before_taxesSS = ((RKSS - 1.0) ./ (1.0 - (TkSS - 1.0))) + 1.0

# resource constaint
CSS = YSS - ISS
