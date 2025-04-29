@R$1 #no replication
# Set aggregate steady state variabel values
A$SS = 1.0
ZI$SS = 1.0
μ$SS = m_par.μ
μw$SS = m_par.μw
Tprog$SS = m_par.Tprog
Tlev$SS = m_par.Tlev
Tc$SS = m_par.Tc

mc$SS = 1.0 ./ μ$SS
mcw$SS = 1.0 ./ μw$SS
Z$SS = m_par.Z
wF$SS = wage(mc$SS, Z$SS, K$SS, N$SS, m_par)
Y$SS = output(Z$SS, K$SS, N$SS, m_par)

σ$SS = 1.0
Tprog_obs$SS = 1.0
Gshock$SS = 1.0
Rshock$SS = 1.0
Tprogshock$SS = 1.0

π$SS = 1.0
πw$SS = 1.0

Sshock$SS = 1.0
RB$SS = m_par.RRB .* π$SS
LP$SS = 1 + RKSS - RBSS
LPXA$SS = 1 + RKSS - RBSS
I$SS = m_par.δ_0 * KSS

BDSS = -sum(distr_bSS .* (n_par.grid_b .< 0) .* n_par.grid_b)

# Calculate taxes and government expenditures
T$SS =
    (Tbar$SS .- 1.0) .* (1.0 ./ m_par.μw .* wF$SS .* N$SS) + # labor income
    (Tbar$SS .- 1.0) .* Π_E$SS + # profit income
    (Tbar$SS .- 1.0) * ((1.0 .- 1.0 ./ m_par.μw) .* wF$SS .* N$SS) # union profit income

RL$SS = RB$SS
RD$SS = RRD$SS .* π$SS

qΠ$SS = m_par.ωΠ .* (1.0 .- 1.0 ./ m_par.μ) .* Y$SS ./ (RRL$SS .- 1 .+ m_par.ιΠ) + 1.0
Bgov$SS = B$SS .- qΠ$SS .+ 1.0
G$SS = T$SS - (RRL$SS - 1.0) * Bgov$SS
Π_F$SS = (1.0 - mc$SS) .* Y$SS
qΠlag$SS = qΠ$SS

C$SS = (Y$SS - m_par.δ_0 * K$SS - G$SS - (RRD$SS .- RRL$SS) * BD$SS)

mcw$SS = 1.0 ./ m_par.μw
u$SS = 1.0

BY$SS = B$SS / Y$SS
TY$SS = T$SS / Y$SS
Tlag$SS = T$SS

Ylag$SS = Y$SS
Bgovlag$SS = Bgov$SS
Glag$SS = G$SS
Ilag$SS = I$SS
wFlag$SS = wF$SS
qlag$SS = q$SS
Clag$SS = C$SS
Tbarlag$SS = Tbar$SS
Tproglag$SS = Tprog$SS

Ygrowth$SS = 1.0
Bgovgrowth$SS = 1.0
Igrowth$SS = 1.0
wgrowth$SS = 1.0
Cgrowth$SS = 1.0
Tgrowth$SS = 1.0

TotalAssets$SS = B$SS + q$SS * K$SS

τprog$SS = Tprog$SS .- 1.0
