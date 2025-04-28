@R$1 #no replication
# Set aggregate steady state variabel values
Tprog$SS = m_par.Tprog
Tlev$SS = m_par.Tlev
Tc$SS = m_par.Tc

mc$SS = 1.0 ./ m_par.μ
mcw$SS = 1.0 ./ m_par.μw
Z$SS = m_par.Z
wF$SS = wage(mc$SS, Z$SS, K$SS, N$SS, m_par)
Y$SS = output(Z$SS, K$SS, N$SS, m_par)

σ$SS = 1.0
Rshock$SS = 1.0

π$SS = 1.0
πw$SS = 1.0

RB$SS = m_par.RRB .* π$SS
I$SS = m_par.δ_0 * KSS

BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)

# Calculate taxes and government expenditures
T$SS =
    (Tbar$SS .- 1.0) .* (1.0 ./ m_par.μw .* wF$SS .* N$SS) + # labor income
    (Tbar$SS .- 1.0) .* Π_E$SS + # profit income
    (Tbar$SS .- 1.0) * ((1.0 .- 1.0 ./ m_par.μw) .* wF$SS .* N$SS) # union profit income

RL$SS = RB$SS
RD$SS = RRD$SS .* π$SS

Bgov$SS = B$SS
G$SS = T$SS - (RRL$SS - 1.0) * Bgov$SS
Π_F$SS = (1.0 - mc$SS) .* Y$SS

C$SS = (Y$SS - m_par.δ_0 * K$SS - G$SS - (RRD$SS .- RRL$SS) * BD$SS)

mcw$SS = 1.0 ./ m_par.μw

Bgovlag$SS = Bgov$SS
wFlag$SS = wF$SS
qlag$SS = q$SS
Tbarlag$SS = Tbar$SS

TotalAssets$SS = B$SS + q$SS * K$SS
