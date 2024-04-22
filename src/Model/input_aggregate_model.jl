#------------------------------------------------------------------------------
# THIS FILE CONTAINS THE "AGGREGATE" MODEL EQUATIONS, I.E. EVERYTHING  BUT THE 
# HOUSEHOLD PLANNING PROBLEM. THE lATTER IS DESCRIBED BY ONE EGM BACKWARD STEP AND 
# ONE FORWARD ITERATION OF THE DISTRIBUTION.
#
# AGGREGATE EQUATIONS TAKE THE FORM 
# F[EQUATION NUMBER] = lhs - rhs
#
# EQUATION NUMBERS ARE GENEREATED AUTOMATICALLY AND STORED IN THE INDEX STRUCT
# FOR THIS THE "CORRESPONDING" VARIABLE NEEDS TO BE IN THE LIST OF STATES 
# OR CONTROLS.
#------------------------------------------------------------------------------

# Magic comment for number (n_rep=2) of code repetitions to cover multiple economies/sectors. 
# $ is the symbol to be repolaced by the number of the economy
# HANK economy is economy 1, rep agent economy is economy 2
# economy 1 does not get a number, economy 2 gets a 2, etc
# entire code gets copied n_rep times with placeholder ($) replaced by the number of the economy
@R$2
#------------------------------------------------------------------------------
# AUXILIARY VARIABLES ARE DEFINED FIRST
#------------------------------------------------------------------------------
# Remaining auxiliary variables

ιΠ = (1.0 ./ 40.0 - 1.0 ./ 800.0) .* m_par.shiftΠ .+ 1.0 ./ 800.0
ωΠ = ιΠ ./ m_par.ιΠ .* m_par.ωΠ

# Elasticities and steepness from target markups for Phillips Curves
η$ = μ$ / (μ$ - 1.0)                                 # demand elasticity
κ$ = η$ * (m_par.κ / m_par.μ) * (m_par.μ - 1.0)     # implied steepness of phillips curve
ηw$ = μw$ / (μw$ - 1.0)                               # demand elasticity wages
κw$ = ηw$ * (m_par.κw / m_par.μw) * (m_par.μw - 1.0) # implied steepness of wage phillips curve

# Capital Utilization
MPKSS = exp(XSS[indexes.rSS]) - 1.0 + m_par.δ_0       # stationary equil. marginal productivity of capital
δ_1 = MPKSS                                        # normailzation of utilization to 1 in stationary equilibrium
δ_2 = δ_1 .* m_par.δ_s                              # express second utilization coefficient in relative terms
# Auxiliary variables
Kserv$ = K$ * u$                                         # Effective capital
MPKserv$ = interest(Kserv$, mc$.*Z$, N$, m_par) .+ m_par.δ_0 # mc .* Z .* m_par.α .* (Kserv ./ N) .^ (m_par.α - 1.0)      # marginal product of Capital
depr$ = m_par.δ_0 + δ_1 * (u$ - 1.0) + δ_2 / 2.0 * (u$ - 1.0)^2.0   # depreciation

Wagesum$ = N$ * w$                                         # Total wages in economy t
Wagesum$Prime = N$Prime * w$Prime                               # Total wages in economy t+1

YREACTION$ = Ygrowth$                                  # Policy reaction function to Y

distr_y = sum(distrSS, dims = (1, 2))

# tax progressivity variabels used to calculate e.g. total taxes
# het agent country
tax_prog_scale = (m_par.γ + m_par.τprog ) / ((m_par.γ + τprog))                        # scaling of labor disutility including tax progressivity
incgross = ((n_par.grid_y ./ n_par.H) .^ tax_prog_scale .* mcw .* w .* N ./ Ht)  # capital liquidation Income (q=1 in steady state)
incgross[end] = (n_par.grid_y[end] .* profits)                         # gross profit income
inc = τlev .* (incgross .^ (1.0 .- τprog))                                 # capital liquidation Income (q=1 in steady state)
taxrev = incgross .- inc                                                 # tax revenues

TaxAux = dot(distr_y, taxrev)
IncAux = dot(distr_y, incgross)

# rep agent country
IncAux2 = mcw2*w2*N2 + profits2
TaxAux2 = (mcw2*w2*N2 + profits2) - τlev2 * (w2*N2 + profits2) .^ (1.0 - τprog2)

Htact = dot(
    distr_y[1:end-1],
    (n_par.grid_y[1:end-1] / n_par.H) .^ ((m_par.γ + m_par.τprog ) / (m_par.γ + τprog)),
)
############################################################################
#           Error term calculations (i.e. model starts here)          #
############################################################################

#-------- States -----------#
# Error Term on exogeneous States
# Shock processes

F[indexes.Gshock$]       = log.(Gshock$Prime) - m_par.ρ_Gshock * log.(Gshock$)     # primary deficit shock
F[indexes.Tprogshock$]   = log.(Tprogshock$Prime) - m_par.ρ_Pshock * log.(Tprogshock$) # tax shock

F[indexes.Rshock$]       = log.(Rshock$Prime) - m_par.ρ_Rshock * log.(Rshock$)     # Taylor rule shock
F[indexes.Sshock$]       = log.(Sshock$Prime) - m_par.ρ_Sshock * log.(Sshock$)     # uncertainty shock

# Stochastic states that can be directly moved (no feedback)
F[indexes.A$]            = log.(A$Prime) - m_par.ρ_A * log.(A$)                # (unobserved) Private bond return fed-funds spread (produces goods out of nothing if negative)
F[indexes.Z$]            = log.(Z$Prime) - m_par.ρ_Z * log.(Z$)                # TFP
F[indexes.ZI$]           = log.(ZI$Prime) - m_par.ρ_ZI * log.(ZI$)             # Investment-good productivity

F[indexes.μ$]            = log.(μ$Prime ./ m_par.μ) - m_par.ρ_μ * log.(μ$ ./ m_par.μ)      # Process for markup target
F[indexes.μw$]           = log.(μw$Prime ./ m_par.μw) - m_par.ρ_μw * log.(μw$ ./ m_par.μw)   # Process for w-markup target

# Endogeneous States (including Lags)
F[indexes.σ] =                                      # only in the het agent economy
    log.(σPrime) -
    (m_par.ρ_s * log.(σ) + (1.0 - m_par.ρ_s) * m_par.Σ_n * log(Ygrowth) + log(Sshock))                     # Idiosyncratic income risk (contemporaneous reaction to business cycle)

F[indexes.Ylag$] = log(Ylag$Prime) - log(Y$)
F[indexes.Bgovlag$] = log(Bgovlag$Prime) - log(Bgov$)
F[indexes.Ilag$] = log(Ilag$Prime) - log(I$)
F[indexes.wlag$] = log(wlag$Prime) - log(w$)
F[indexes.Tlag$] = log(Tlag$Prime) - log(T$)
F[indexes.qlag$] = log(qlag$Prime) - log(q$)
F[indexes.Clag$] = log(Clag$Prime) - log(C$)
F[indexes.av_tax_ratelag$] = log(av_tax_ratelag$Prime) - log(av_tax_rate$)
F[indexes.τproglag$] = log(τproglag$Prime) - log(τprog$)
F[indexes.qΠlag$] = log(qΠlag$Prime) - log(qΠ$)
F[indexes.plag$] = log(plag$Prime) - log(p$)

# Growth rates
F[indexes.Ygrowth$] = log(Ygrowth$) - log(Y$ / Ylag$)
F[indexes.Tgrowth$] = log(Tgrowth$) - log(T$ / Tlag$)
F[indexes.Bgovgrowth$] = log(Bgovgrowth$) - log(Bgov$ / Bgovlag$)
F[indexes.Igrowth$] = log(Igrowth$) - log(I$ / Ilag$)
F[indexes.wgrowth$] = log(wgrowth$) - log(w$ / wlag$)
F[indexes.Cgrowth$] = log(Cgrowth$) - log(C$ / Clag$)

#  Taylor rule (PPI based) and interest rates
F[indexes.RB$] =
    log(RB$Prime) - XSS[indexes.RB$SS] - ((1 - m_par.ρ_R) * m_par.θ_π) .* log(π$) -
    ((1 - m_par.ρ_R) * m_par.θ_Y) .* log(YREACTION$) -
    m_par.ρ_R * (log.(RB$) - XSS[indexes.RB$SS]) - log(Rshock$)

# Tax rule
F[indexes.τprog$] =
    log(τprog$) - m_par.ρ_P * log(τproglag$) - (1.0 - m_par.ρ_P) * (XSS[indexes.τprog$SS]) -
    (1.0 - m_par.ρ_P) * m_par.γ_YP * log(YREACTION$) -
    (1.0 - m_par.ρ_P) * m_par.γ_BP * (log(Bgov$) - XSS[indexes.Bgov$SS]) - log(Tprogshock$)


F[indexes.τlev$] = av_tax_rate$ - TaxAux$ ./ IncAux$  # Union profits are taxed at average tax rate
F[indexes.T$] = log(T$) - log(TaxAux$ + av_tax_rate$ * unionprofits$)


F[indexes.av_tax_rate$] =
    log(av_tax_rate$) - m_par.ρ_τ * log(av_tax_ratelag$) -
    (1.0 - m_par.ρ_τ) * XSS[indexes.av_tax_rate$SS] -
    (1.0 - m_par.ρ_τ) * m_par.γ_Yτ * log(YREACTION$) -
    (1.0 - m_par.ρ_τ) * m_par.γ_Bτ * (log(Bgov$) - log(Bgovlag$))#XSS[indexes.BgovSS])

# --------- Controls ------------
# Deficit rule
F[indexes.π$] =
    log(Bgovgrowth$Prime) + m_par.γ_B * (log(Bgov$) - XSS[indexes.Bgov$SS]) -
    m_par.γ_Y * log(YREACTION$) - m_par.γ_π * log(π$) - log(Gshock$)

F[indexes.G$] = log(p$*G$) - log(Bgov$Prime + T$ - RB$ / πCPI$ * Bgov$)             # Government Budget Constraint, perfect home bias of G

# Phillips Curve to determine equilibrium markup, output, factor incomes 
F[indexes.mc$] =
    (log.(π$) - XSS[indexes.π$SS]) - κ$ * (mc$ - 1 ./ μ$) -
    m_par.β * ((log.(π$Prime) - XSS[indexes.π$SS]) .* Y$Prime ./ Y$)

# Wage Phillips Curve 
F[indexes.mcw$] =
    (log.(πw$) - XSS[indexes.πw$SS]) - (
        κw$ * (mcw$ - 1 ./ μw$) +
        m_par.β * ((log.(πw$Prime) - XSS[indexes.πw$SS]) .* Wagesum$Prime ./ Wagesum$)
    )
# worker's wage = mcw * firm's wage

# Wage Dynamics
F[indexes.πw$] = log.(w$ ./ wlag$) - log.(πw$ ./ πCPI$)                   # Definition of real wage inflation

# Capital utilisation
F[indexes.u$] = MPKserv$ - q$ * (δ_1 + δ_2 * (u$ - 1.0))           # Optimality condition for utilization

# Prices
F[indexes.r$] = log.(r$) - log.(1 + MPKserv$ * u$ - q$ * depr$)       # rate of return on capital

F[indexes.mcww$] = log.(mcww$) - log.(mcw$ * w$)                        # wages that workers receive

F[indexes.w$] = log.(w$) - log.(wage(Kserv$, Z$ * mc$, N$, m_par))     # wages that firms pay

F[indexes.unionprofits$] = log.(unionprofits$) - log.(w$ .* N$ .* (1.0 - mcw$))  # profits of the monopolistic unions

# firm_profits: price setting profits + investment profits. The latter are zero and do not show up up to first order (K'-(1-δ)K = I).
F[indexes.firm_profits$] =
    log.(firm_profits$) - log.(Y$ .* (p$ - mc$) .+ q$ .* (K$Prime .- (1.0 .- depr$) .* K$) .- I$)
F[indexes.profits$] = log.(profits$) - log.((1.0 .- ωΠ) .* firm_profits$ .+ ιΠ .* (qΠ$ .- 1.0)) # distributed profits to entrepreneurs
F[indexes.qΠ$] =
    log.(RB$Prime ./ πCPI$Prime) .-
    log.(((qΠ$Prime .- 1.0) .* (1 - ιΠ) .+ ωΠ .* firm_profits$Prime) ./ (qΠ$ .- 1.0))

F[indexes.RL] = log.(RL) - log.((RB .* Bgov .+ 
                    πCPI .* ((qΠ .- 1.0) .* (1 - ιΠ) .+ ωΠ .* firm_profits) .+
                    B12./rer .* RB2 .* πCPI ./ πCPI2) ./ B)

F[indexes.RL2] = log.(RL2) - log.((RB2 .* Bgov2 .+ 
                    πCPI2 .* ((qΠ2 .- 1.0) .* (1 - ιΠ) .+ ωΠ .* firm_profits2) .-
                    RB2 .* B12 .* (m_par.α_S / (1.0 - m_par.α_S))) ./ B2)

F[indexes.Bgov] = log.(B) - log.(Bgov + (qΠlag .- 1.0) + B12./rerlag)                                 # total liquidity demand
F[indexes.Bgov2] = log.(B2) - log.(Bgov2 + (qΠlag2 .- 1.0) - B12 .* (m_par.α_S / (1.0 - m_par.α_S)))  

F[indexes.q$] =
    1.0 - ZI$ * q$ *
    (1.0 - m_par.ϕ / 2.0 * (Igrowth$ - 1.0)^2.0 - # price of capital investment adjustment costs
     m_par.ϕ * (Igrowth$ - 1.0) * Igrowth$) -
    m_par.β * ZI$Prime * q$Prime * m_par.ϕ * (Igrowth$Prime - 1.0) * (Igrowth$Prime)^2.0

# Asset market premia
F[indexes.LP$] = log.(LP$) - (log((q$ + r$ - 1.0) / qlag$) - log(RL$ / πCPI$))                   # Ex-post liquidity premium           
F[indexes.LPXA$] = log.(LPXA$) - (log((q$Prime + r$Prime - 1.0) / q$) - log(RL$Prime / πCPI$Prime))  # ex-ante liquidity premium

# Aggregate Quantities
F[indexes.I$] =
    K$Prime .- K$ .* (1.0 .- depr$) .-
    ZI$ .* I$ .* (1.0 .- m_par.ϕ ./ 2.0 .* (Igrowth$ - 1.0) .^ 2.0)           # Capital accumulation equation

F[indexes.N$] = log.(N$) - log.(labor_supply(w$.*mcw$, τlev$, τprog$, Ht$, m_par)) # Labor supply equation

F[indexes.Y$] = log.(Y$) - log.(output(Kserv$, Z$, N$, m_par))            # Output equation

# Error Term on prices/aggregate summary vars (logarithmic, controls), here difference to SS value averages
F[indexes.BY$] = log.(BY$) - log.(B$ / Y$)                                                               # Bond to Output ratio
F[indexes.TY$] = log.(TY$) - log.(T$ / Y$)                                                               # Tax to output ratio

# Open Economy equations

## International prices

F[indexes.p12]      = log(p2)  - log(p12 * rer) 
F[indexes.p21]      = log(p21) - log(p * rer) 

# Dynamic LOOP
F[indexes.rer]   =  log(rer/der)  - log(πCPI2 / πCPI *  p12lag/ p2lag) 
    
F[indexes.der]   =  log(RBPrime/RB2Prime) - log(derPrime) #der: delta exchange rate
#'CPI index'
F[indexes.p]   = log(1.0)       - log(((1.0 - (1.0 - m_par.α_S) * m_par.ω) * p^(1.0 - m_par.ϵ_e) + (1.0 - m_par.α_S) * m_par.ω * p12^(1 - m_par.ϵ_e))^(1.0/(1.0 - m_par.ϵ_e)))
F[indexes.p2]   = log(1.0)       - log((m_par.α_S * m_par.ω * p21^(1.0 - m_par.ϵ_e) + (1.0 - m_par.α_S * m_par.ω) * p2^(1.0 - m_par.ϵ_e))^(1.0/(1.0 - m_par.ϵ_e)))

#'PPI Inflation Definition
F[indexes.πCPI$]  =  log(p$/plag$)      - log((π$/πCPI$))               


#'Home net exports'
F[indexes.nx]          =  p * (Y - G - nx)  - (C + I + m_par.Rbar*BD - (A - 1.0) * (RL * B / πCPI))
F[indexes.C]           =  log(Y)   - log(p^(-m_par.ϵ_e) * ((1 - (1 - m_par.α_S) * m_par.ω) * 
                                            (C + I + BD*m_par.Rbar - (A - 1.0) * (RL * B / πCPI)) + 
                                            (1 - m_par.α_S) * m_par.ω * rer^(-m_par.ϵ_e) * (C2 + I2 + BD2*m_par.Rbar - 
                                            (A2 - 1.0) * (RL2 / πCPI2 * B2 ))) + 
                                            G)
F[indexes.C2]            =  log(Y2)   - log(p2^(-m_par.ϵ_e) * (m_par.α_S * m_par.ω * rer^(m_par.ϵ_e) * 
                                            (C + I+ BD*m_par.Rbar - (A - 1.0) * (RL * B / πCPI)) +
                                            (1 - m_par.α_S * m_par.ω) * (C2 + I2 + BD2*m_par.Rbar - 
                                            (A2 - 1.0) * (RL2 / πCPI2 * B2))) + 
                                            G2)


# Distribution summary statistics used in this file (using the steady state distrubtion in case). 
# Lines here generate a unit derivative (distributional summaries do not change with other aggregate vars).

# Heterogeneous agent economy
F[indexes.K]  = log.(K) - XSS[indexes.KSS]                                                       # Capital market clearing           
F[indexes.B]  = log.(B) - XSS[indexes.BSS]                                                       # Bond market clearing
F[indexes.BD] = log.(BD) - XSS[indexes.BDSS]                                                     # IOUs

F[indexes.Ht] = log.(Ht) - log.(Htact)
F[indexes.τlev] = av_tax_rate - TaxAux ./ IncAux
F[indexes.T] = log(T) - log(TaxAux + av_tax_rate * unionprofits)

# Representative agent economy (economy 2 ...)
F[indexes.K2] = log(LP2) - log(m_par.LP2)                                                       # fixed spread on capital returns and bond returns in 2
F[indexes.B2] = log((C2 - 1.0./(1+m_par.γ)*N2.^(1.0 .+ m_par.γ)).^(-m_par.ξ)) -                 # consumption Euler equation for country 2
                m_par.β.*RLPrime2/πCPI2Prime.*log((C2Prime - 1.0./(1+m_par.γ)*N2Prime.^(1.0 .+ m_par.γ)).^(-m_par.ξ))

F[indexes.B12] = log((1.0 - av_tax_rate) * (w*N + profits) + #all household incomes are consumed or saved, B12 savings odf homae abroad (in liquid assets)
                     (p* Y - w*N - profits) + Tr + A * B * RL/πCPI)    - 
                 log(C + I + BD * m_par.Rbar  + BgovPrime + (qΠ .- 1.0) + B12Prime/rer) 

# Add distributional summary stats that do change with other aggregate controls/prices and with estimated parameters
F[indexes.Ht2] = log.(Ht2) - log.(1.0)
F[indexes.T2]  = log(T2) - log(TaxAux2 + av_tax_rate2 * unionprofits2)
F[indexes.τlev2] = av_tax_rate2 - TaxAux2 ./ IncAux2

# other distributional statistics not used in other aggregate equations and not changing with parameters, 
# but potentially with other aggregate variables are NOT included here. They are found in FSYS.
