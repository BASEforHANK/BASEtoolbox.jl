@R$1 #no replication

Z$SS = m_par.Z

mcw = 1.0 ./ m_par.μw
mc = 1.0 ./ m_par.μ

wF$SS = wage(mc, Z$SS, K$SS, N$SS, m_par)
Y$SS = output(Z$SS, K$SS, N$SS, m_par)
Tprog$SS = m_par.Tprog
Tlev$SS = m_par.Tlev

I$SS = m_par.δ_0 * KSS
σSS = m_par.σ
qSS = 1.0
Tc$SS = m_par.Tc
BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)

# Calculate taxes and government expenditures
RRD$SS = RRL$SS

Π_F$SS = (1.0 - mc) .* Y$SS

C$SS = Y$SS - I$SS

TotalAssets$SS = K$SS
