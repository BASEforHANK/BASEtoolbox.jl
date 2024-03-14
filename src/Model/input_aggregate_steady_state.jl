@R$1 #no replication 
# Set aggregate steady state variabel values
A$SS = 1.0
Z$SS = 1.0
ZI$SS = 1.0
μ$SS = m_par.μ
μw$SS = m_par.μw
τprog$SS = m_par.τprog 
τlev$SS = m_par.τlev

σ$SS = 1.0
τprog_obs$SS = 1.0
Gshock$SS = 1.0
Rshock$SS = 1.0
Tprogshock$SS = 1.0

Sshock$SS = 1.0
RB$SS = m_par.RB
LP$SS = 1 + rSS - RBSS
LPX$ASS = 1 + rSS - RBSS
I$SS = m_par.δ_0 * KSS

π$SS = 1.0
πw$SS = 1.0

BDSS = -sum(distr_mSS .* (n_par.grid_m .< 0) .* n_par.grid_m)
# Calculate taxes and government expenditures
T$SS = dot(distr_y$SS, taxrev) + av_tax_rate$SS * ((1.0 .- 1.0 ./ m_par.μw) .* w$SS .* N$SS)

Bgov$SS = B$SS .- qΠSS_fnc(Y$SS, m_par.RB, m_par) .+ 1.0
G$SS = T$SS - (m_par.RB ./ m_par.π - 1.0) * Bgov$SS
mc$SS = 1.0 ./ m_par.μ
firm_profits$SS = (1.0 - mc$SS) .* Y$SS
qΠ$SS = qΠSS_fnc(Y$SS, RB$SS, m_par)
qΠlag$SS = qΠ$SS
RL$SS = m_par.RB

C$SS = (Y$SS - m_par.δ_0 * K$SS - G$SS - m_par.Rbar * BD$SS)

q$SS = 1.0
mcw$SS = 1.0 ./ m_par.μw
mcww$SS = w$SS * mcw$SS
u$SS = 1.0
unionprofits$SS = (1.0 - mcw$SS) .* w$SS .* N$SS

BY$SS = B$SS / Y$SS
TY$SS = T$SS / Y$SS
Tlag$SS = T$SS

Ylag$SS = Y$SS
Bgovlag$SS = Bgov$SS
Glag$SS = G$SS
Ilag$SS = I$SS
wlag$SS = w$SS
qlag$SS = q$SS
Clag$SS = C$SS
av_tax_ratelag$SS = av_tax_rate$SS
τproglag$SS = τprog$SS

Ygrowth$SS = 1.0
Bgovgrowth$SS = 1.0
Igrowth$SS = 1.0
wgrowth$SS = 1.0
Cgrowth$SS = 1.0
Tgrowth$SS = 1.0
Ht$SS = 1.0

# calculate and display aggregate moments
println("B$SS/Y$SS: ", B$SS / Y$SS)
println("Stockshare$SS: ", (qΠSS_fnc(Y$SS, m_par.RB, m_par) .- 1.0) / B$SS)
println("Bgov$SS/Y$SS: ", Bgov$SS / Y$SS)
println("G$SS/Y$SS: ", G$SS / Y$SS)
println("firm_profits$SS: ", firm_profits$SS)
