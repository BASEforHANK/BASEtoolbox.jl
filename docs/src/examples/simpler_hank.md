# Simpler HANK model
This model is a simplified version of the model in [Bayer, Born, and Luetticke (2024, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20201875), henceforth BBL (2024), which is implemented as a detailed [baseline example](#baseline.md).
It combines a downsized version of a business cycle model in the style of Smets and Wouters (2007) with an incomplete markets model of heterogeneous agents and portfolio choice as in [Bayer et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA13601).

## Economics

The general equilibrium model includes a firm sector, government sector and household sector composed of heterogeneous agents, who make portfolio decision in incomplete markets. The two assets in the household portfolio are liquid bonds and capital, which is held as illiquid asset as in [Bayer et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA13601). Households supply labor under idiosyncratic labor income risk. This idiosyncratic labor income risk is time-invariant in the sense that we keep the variance of idiosyncratic shocks fixed over time in this simple model version. Raw labor is differentiated by labor unions, allowing the modelling of Calvo-style nominal wage stickiness with a wage Phillip’s curve. The union's labor services are further bundeled into labor services by a competitive labor packer. While final goods are produced competitively, intermediate goods are produced under monopolistic competition and Calvo price rigidity. The arising firm profits are consumed by entrepreneus only, whereas households are either entrepreneurs or workers. Workers in turn receive the profits from labor unions. The standard New Keynesian public sector (with a fiscal authority and a central bank), is extended by progressive taxation on labor incomes and profits, with a fixed form of progressivity.
To simplify, this simple version of the model does not include flexible capital utilization or capital adjustment costs. Moreover, the model features a very reduced set of shocks compared to standard DSGE business cycle models. The only shocks included in this simplified model are TFP technology shocks and monetary policy shocks.

### Details

#### Households

Households transition stochastically between the roles of entrepreneurs and workers and face idiosyncratic labor productivity risk in this simplified version of the BBL (2024) model. Entrepreneurs receive and consume all firm profits, while workers earn labor income and union profits. Given their incomes households optimize their consumption and savings intertemporally, subject to a budget constraint. They allocate their savings between capital as illiquid asset, and liquid assets. Capital is considered illiquid because households can only participate in the capital market with a certain probability each period. Households that do not participate in the capital market ($k_{i,t+1} = k_{i,t}$) still receive dividends and can adjust their liquid asset holdings. Borrowing in liquid assets is constrained by an exogenous debt limit $\underline{B}$, and the return on liquid assets depends on whether the household is a borrower or a lender.
Please see the documentation of the [household problem](../HouseholdProblem.md) for further details and all household equations deeper in the package's code.

#### Government Sector


The government sector entails a monetary and a fiscal authority.

##### Fiscal Policy

There is tax progressivity on labor income and profits, such that individual after-tax income is:

```math
\begin{align*}
y_{it} = (1- \tau^L_t)(w_t h_{it} n_{it})^{1- \tau^P}
\end{align*}
```
Relative to BBL (2024) we keep the degree of progressivity $\tau^P$ fixed over time.

The fiscal authority issues government bonds and levies taxes according to two separate rules. The first rule defines bond issuing and is similar to equation (33) in BBL (2024) without structural deficit shocks:
```math
\begin{align*}
\frac{B^{gov}_{t+1}}{B^{gov}_{t}} = \left( \frac{B^{gov}_t}{\bar{B}^{gov}}\right)^{-\gamma_B} \left( \frac{\pi_t}{\bar{\pi}} \right)^{\gamma_{\pi}} \left( \frac{Y_t}{\bar{Y}} \right)^{\gamma_Y}. \tag{indexes.$\pi$}
\end{align*}
```

The second rule pins down the average tax rate $\tau_t$ in the economy and is equivalent to equation (34) in BBL (2024):
```math
\begin{align*}
\frac{\tau_t}{\bar{\tau}} = \left( \frac{\tau_{t-1}}{\bar{\tau}}\right)^{\rho_\tau}  \left( \frac{B^{gov}_t}{B^{gov}_{t-1}} \right)^{(1-\rho_\tau) \gamma_B^\tau} \left( \frac{Y_t}{\bar{Y}} \right)^{(1-\rho_\tau) \gamma_Y^\tau}. \tag{indexes.Tbar}
\end{align*}
```

The tax level parameter $\tau^L_t$ is then determined such that the average tax rate $\tau_t$ is obtained under the progressive taxing system, following equation (35) in BBL (2024):
```math
\begin{align*}
\tau_t =
\tau_t \left(
    \ \frac{w^H_t \cdot N_t}{H_{prog,t}},\ \tau^L_{t},\ \tau^p,\ \Pi^E_t,\ \tilde{H}_{t},\ \text{distr}^h_t \tag{indexes.Tlev}
\right)
\end{align*}
```
!!! note
    Here we can rely on a function used in the household problem, $\tau_t(\cdot)$ that computes the average tax rate given the variables for the income taxation, $\tau^L_t$ and $\tau^p$ (constant here), labor income, and entrepreneurs' profits as well as the distribution about productivity types.

Total taxes can then be calculated based on the average tax rate:
```math
\begin{align*}
T_t = \tau_t w^H_t N_t + \tau_t \Pi^E_t + \tau_t \Pi^U_t \tag{indexes.T}
\end{align*}
```

Government spending is given residually through the government's budget constraint:
```math
\begin{align*}
G_t = B^{gov}_{t+1} + T_t - \frac{R^B_t}{\pi_t} B^{gov}_{t} \tag{indexes.G}
\end{align*}
```

##### Monetary Policy

The monetary authority follows a Taylor Rule when setting the nominal interest rate $RB_t$. It stabilizes inflation and output growth while smoothing interest rates over time similar to equation (32) in BBL (2024):

```math
\begin{align*}
\frac{R^B_{t+1}}{\bar{R}^B} =  \pi_t ^{(1-\rho_R) \theta_\pi} \left(\frac{Y_t}{\bar{Y}} \right) ^{(1-\rho_R) \theta_Y} \left( \frac{R^B_t}{\bar{R}^B} \right) ^ {\rho_R} \epsilon^R_t \tag{indexes.RB}
\end{align*}
```
The monetary policy shock $\epsilon^R_t$ follows a standard AR(1)-process:
```math
\begin{align*}
\log \epsilon^R_{t+1} = \rho_R \log \epsilon^R_{t} + \tilde{\epsilon}^R_{t+1} \tag{indexes.Rshock}
\end{align*}
```


#### Firm Sector

The firm sector of this simple example closely follows the structure of the business cycle models as in Smets and Wouters (2007), but under a number of simplifications.
Importantly, this simple version of the model does not include flexible capital utilization or capital adjustment costs.

As described in the beginning of section I.A of BBL (2024), there are final goods producer, who bundle differentiated intermediate goods with a Dixit-stiglitz aggregator and elasticity of substitution $\eta_t$.

```math
Y_t = \left( \int y_{jt}^{\frac{\eta - 1}{\eta}} \, dj \right)^{\frac{\eta}{\eta - 1}},
```
Given price $p_{jt}$ of an intermediate good the demand for it is:

```math
y_{jt} = \left( \frac{p_{jt}}{P_t} \right)^{-\eta} Y_t

```

with $P_t$ being given by $P_t = \left( \int p_{jt}^{1 - \eta} \, dj \right)^{1 / (1 - \eta)}$.

In this simplified version of BBL (2024) there is no variable capital utilization and the production function of intermediate goods producers is:
```math
\begin{align*}
y_{jt} = Z_t N_{jt}^{1-\alpha} K_{jt}^{\alpha}
\end{align*}
```

Consequently, there is linear depreciation (instead of depreciation dependant on capital utilization) and $\delta_0$ denotes the constant depreciation rate of capital goods.
Moreover, total factor productivity $Z_t$ follows an AR(1) process in logs.
Intermediate goods producers minimize costs to meet their demand
```math
\begin{align*}
w^F_t N_t + (r_t + q_t δ_0) K_t,
\end{align*}
```
where the last term is rental rate of capital plus the cost of capital depreciation $(r_t + q_t δ_0)$ and $w^F_t$ is the real wage paid to the labor packers.
With factor markets being competitive, the first order conditions are:
```math
\begin{align*}
w_t^F &= (1-\alpha) \, mc_{t} Z_t \left( \frac{K_t}{N_t} \right)^{\alpha} \\
&= \text{wage} \left( mc_t, \ Z_t, \ K_t, \ N_t \right)
\tag{indexes.wF} \\
r_t + q_t \delta_0
&= \alpha mc_{t} Z_t \left( \frac{N_{t}}{ K_{t}} \right)^{1-\alpha}.
\end{align*}
```
Marginal costs $mc_t$ are constant across firms due to symmetry and the constant returns to scale production function.

Intermediate goods producers face price adjustment frictions à la Calvo (1983) and maximize the present value of real profits, as in BBL (2024). Prices are indexed to steady state inflation and $1-\lambda_Y$ describes the probability of adjustment. This results in the firm's Phillip's Curve of the gross inflation rate $\pi_t:=P_t/P_{t-1}$:
```math
\begin{align*}
    \log \pi_{t} - \log \bar{\pi} = \kappa \left(mc_t - \frac{1}{\mu} \right) + \beta \left(( \log \pi_{t+1} - \log \bar{\pi}) \frac{Y_{t+1}}{Y_t} \right)
,
\tag{indexes.mc}
\end{align*}
```
with all irrelevant terms for a first-order approximation dropped and $\kappa = \frac{(1 - \lambda_Y)(1 - \lambda_Y \beta)}{\lambda_Y}$.

 We assume the firm's time constant discount factor $\beta$ to be equal to the discount factor of households, see BBL (2024), p. 1216, for a setup with risk-neutral managers justifying this.
 In this simplified model, the target markups are constant both for the production price and wages ($\mu$ and $\mu_w$ are constant), because the elasiticities of substitution ($\eta$ and $\zeta$) are constant.

 Without capital adjustment costs and constant marginal efficiency of investment the (real) price of capital simplifies to a constant which is 1 (see below).

 The law of motion for aggregate capital is simply:
```math
\begin{align*}
    K_{t+1} = I_t + (1 - \delta_0) K_t \tag{indexes.I}
\end{align*}
```

Real profits of the firms consisting of monopolistic profits by the intermediate goods producers and the profits by producers of capital are as follows:
```math
\begin{align*}
 \Pi^F_t = Y_t(1- mc_t) + q_t (K_{t+1} - (1 - \delta_0) K_t) - I_t
\tag{indexes.$\Pi$\_F}
\end{align*}
```

Given the unitary price of capital and the law of motion for capital, average profits reduce to monopolistic profits alone, which in real terms are:

```math
\begin{align*}
 \Pi^F_t = \frac{Y_t (P_t-MC_t)}{P_t} = Y_t(1- mc_t)
\end{align*}
```

In this example all profits by firms are distributed to the entrepreneurs (i.e. they do not sell any claims to their profits as in BBL (2024)):
```math
\begin{align*}
 \Pi^E_t = \Pi^F_t \tag{indexes.$\Pi$\_E}
\end{align*}
```

The net return on capital $RK_t$ is equal to the rental factor of capital $r_t + 1$, and thus:
```math
\begin{align*}
RK_t &= 1 + r_t = 1 + (1 - \alpha) mc_{t} Z_t \left( \frac{N_{jt}}{ K_{jt}} \right)^{1-\alpha} - q_t \delta_0 \\
&= 1 + \text{interest} \left( mc_t, \ Z_t, \ K_t, \ N_t\right) + \delta_0 - q_t \delta_0.
\tag{indexes.RK}
\end{align*}
```

With symmetric firms, aggregate output can be computed based on aggregates:
```math
\begin{align*}
Y_t &= Z_t K_t^{\alpha} N_t^{1 - \alpha} \\
&= \text{output} \left( Z_t, \ K_t, \ N_t \right).
\tag{indexes.Y}
\end{align*}
```

Aggregate productivity evolves according to an AR(1) process:
```math
\begin{align*}
\log Z_{t+1} = \rho_Z \log Z_{t} + \tilde{\epsilon}^Z_{t+1},
\tag{indexes.Z}
\end{align*}
```
where $\tilde{\epsilon}^Z_{t+1}$ is a TFP shock.

#### Labor Unions and Packers

There is a continuum of labor unions, who differentiate the raw labor supply from households.
Above, we have mentioned that factor markets are competitive. Behind this lies a competitive labor packer, who bundles the differentiated labor services from unions for the intermediate goods producers with a Dixit-Stiglitz aggregator and elasticity of substitution $\zeta$.

```math
N_t = \left( \int \hat{n}_{jt}^{\frac{\zeta - 1}{\zeta}} \, dj \right)^{\frac{\zeta}{\zeta - 1}},
```

The labor packer receives the nominal wage $W_t^F$ for the bundled labor service.
Given a union's nominal wage $W_{jt}$ for its differentiated labor variety, the labor packer's optimal demand for it is:

```math
\hat{n}_{jt} = \left( \frac{W_{jt}}{W_t^F} \right)^{-\zeta} N_t,
```

The unions face price rigidities à la Calvo for their specific wages $W_{jt}$, which gives rise to a New Keynesian Wage Phillips Curve with constant target markup $\mu_w$:
```math
\begin{align*}
    \log \pi^w_t - \bar{\pi}_w = \kappa_w \left(mc^w_t - \frac{1}{\mu_w} \right) + \beta \left( \log \pi^w_{t+1} - \bar{\pi}^w \frac{N_{t+1} w^F_{t+1}}{N_t w^F_t} \right) \tag{indexes.mcw}
\end{align*}
```

where inflation of nominal wages $\pi^w_t$ is defined by:
```math
\begin{align*}
    \pi^w_t = \pi_t \frac{w^F_t}{w^F_{t-1}} \tag{indexes.$\pi$w}
\end{align*}
```

Wages households receive are given by the markdown on the union's received wages:
```math
\begin{align*}
    w^H_t = mc^w_t \cdot w^F_t \tag{indexes.wH}
\end{align*}
```

Average profits of the unions, which are allocated back to the workers, are in real terms:

```math
\begin{align*}
 \Pi^U_t = (w_t^F - w_t^H) N_t \tag{indexes.$\Pi$\_U}
\end{align*}
```

When defining $mc_{w,t}$ as the realized markdown on union's received wages for households, this is equal to
```math
\begin{align*}
    \Pi^U_t = w_t^F N_t (1 - \text{mc}_{w,t} )
\end{align*}
```

Average Labor Productivity
```math
\begin{align*}
    H_{prog, t} =
    \sum_{i=1}^{n-1}
        \text{distr}^h_{i,t} \cdot
        \left( \frac{h_i}{\tilde{H}_t} \right)^{\text{scale}(\tau^p)}
\tag{indexes.Hprog}
\end{align*}
```

#### Asset markets

The nominal lending rate for IOUs has to be equal to the nominal rate set by the monetary authority:
```math
\begin{align*}
    R^L_t = R^B_t \tag{indexes.RL}
\end{align*}
```

Real rates for debt, $r^D_t$, and lending, $r^L_t$, are defined as:
```math
\begin{align*}
    r^D_t &= \frac{R^D_t}{\pi_t} \tag{indexes.RD} \\
    r^L_t &= \frac{R^L_t}{\pi_t} \tag{indexes.RRL}
\end{align*}
```

There is a fixed interest rate spread $\bar{R}$ between the real borrowing rate $r^D_t$ and real lending rate $r^L_t$:
```math
\begin{align*}
    r^D_t &= r^L_t + \bar{R} \\
    &= \text{borrowing\_rate\_ss} \left( r^L_t \right)
    \tag{indexes.RD}
\end{align*}
```

#### Market Clearing

The only liquid asset in the economy is government bonds (different from BBL (2024)):
```math
\begin{align*}
    B_t = B^{gov}_t \tag{indexes.B},
\end{align*}
```
IOUs are in zero net supply.

The final goods market clearing condition or equivalently the aggregate resource constraint is:
```math
\begin{align}
Y_t = G_t + I_t + C_t + B^D_t (r^D_t - r^L_t)
\tag{indexes.C}
\end{align}
```

The total demanded assets of households $B^D_t$ (sum of liquid asset holdings $< 0$) enter this equation because the interest rate spread paid on them is neither available for consumption nor for investment or government expenditure.

Aggregate labor demand must equal labor supply by households:
```math
\begin{align*}
N_t = N^{supply}_t \left( \ w^H_t, \ H_{prog,t}, \ \tau^L_t, \ τ^p, \ τ^C, \ \text{taxbase} \right) \tag{indexes.N}
\end{align*}
```
where we use the function $N^{supply}_t(\cdot)$ from the household problem and insert the $\text{taxbase} = \frac{w^H_t \cdot N_t}{H_{prog,t}} + \Pi^E_t$. Again, in this simple version, the tax rates $\tau^p$ and $\tau^C$ are constant.

We also track the number of total assets by the following accounting identity:
```math
\begin{align*}
\text{total assets} = q_{t-1} K_t + B_t \tag{indexes.TotalAssets}
\end{align*}
```

## Additional equations

### Constants entering the household problem

We set some variables entering the household problem to be constant. In this example, this applies to $\tau^p$, $\tau^C$, $\sigma$, and $q$. This is done by:
```math
\begin{align*}
\tau^p_t &= \bar{\tau}^p \tag{indexes.Tprog} \\
\tau^C_t &= \bar{\tau} \tag{indexes.Tc} \\
\sigma_t &= \bar{\sigma} \tag{indexes.$\sigma$} \\
q &= \bar{q}  \tag{indexes.q}
\end{align*}
```

!!! danger
    The model expects these variables to be defined. Therefore, you cannot simply remove them from the model equations or alter their names.

### Lags

We need to define lags by their variable primes and current values, as explained in the section *Lags and leads* of the [general example strucuture guideline](../GeneralStructure.md). This applies to the variables $B^{gov}_t$, $w^F_t$, $q_t$, and $\tau_t$. Those equations are (indexes.Bgovlag), (indexes.wFlag), (indexes.qlag), and (indexes.Tbarlag).

### Closing the aggregate model

The equations that we add in order to close the aggregate model are the market clearing conditions for capital, bonds and private debt (IOUs), as well as the equation for the average labor productivity, i.e. equations (indexes.K), (indexes.B), (indexes.BD) and (indexes.Htilde).

### Auxiliary variables

In this example we use $\text{distr}_h$ as auxiliary variable, which sums up the distribution along the bond and capital grid to get the number of households for each productivity state.

## Steady state
The steady state of this model economy is characterized by multiple conditions, supplied in an extra `.mod` file.

#### Firm Sector

From the price Phillips Curve and the wage Phillips Curve, it follows that in the steady state, marginal costs are equal to the reciprocal of the target markup, and realized markdowns on wages paid to households by unions are the reciprocal of the wage markup:
```math
\begin{align*}
\bar{mc} &= 1/\mu \\
\bar{mc}^w &= 1/\mu_w.
\end{align*}
```

Price and wage inflation are normalized to one in the steady state. In addition, total factor productivity $Z$ is set to one in the steady state consistent with its AR(1) process in logs.

Output and wages in steady state are:
```math
\begin{align*}
\bar{Y} &= \bar{Z} \bar{K}^\alpha \bar{N}^{1 - \alpha}\\
\bar{w}^F &= \bar{mc} (1 - \alpha) \, \bar{Z} \left( \frac{\bar{K}}{\bar{N}} \right)^{\alpha}
\end{align*}
```

Firm profits in steady state are equal to the monopolisitic profits
```math
\begin{align*}
\bar{\Pi}^F = (1 - \bar{mc}) \bar{Y}
\end{align*}
```

Consistent with the law of motion for aggregate capital, aggregate investment is equal to depreciation in steady state:
```math
\begin{align*}
\bar{I} = \delta_0 \bar{K}
\end{align*}
```

#### Government Sector

The monetary policy shock $\epsilon^R_t$ is set to one in the steady state, consistent with its AR(1) process in logs.

Total tax revenue of the government is composed of labor income, profit income, and union profit income in the steady state:
```math
\begin{align*}
\bar{T} = \overline{\tau} \left( \frac{1}{\mu_w} \overline{w}_F \overline{N} \right)
+ \overline{\tau} \overline{\Pi}_E
+ \overline{\tau}  \left( \left(1 - \frac{1}{\mu_w} \right) \overline{w}_F \overline{N} \right) = \overline{\tau} \overline{\Pi}_E + \overline{\tau} \overline{w}_F \overline{N}
\end{align*}
```

Government spending is given residually in the steady state:
```math
\begin{align*}
\bar{G} = \bar{T} - (\bar{r}^L - 1.0) * \bar{B}^{gov}
\end{align*}
```

#### Asset markets
Steady-state nominal interest rates are calculated from gross real rates (equal to one here) and steady-state inflation:
```math
\begin{align*}
\bar{R}^B &= \bar{r}^B_t \bar{\pi} = \bar{r}^L_t \bar{\pi} = \bar{R}^L \\
\bar{R}^D &= \bar{r}^D_t \bar{\pi}
\end{align*}
```

In the steady state, the total supply of liquid assets is also simply equal to government bonds, because IOUs are in zero net supply.
```math
\begin{align*}
    \bar{B} = \bar{B}^{gov}
\end{align*}
```

Total assets in steady state are given by bonds and the value of capital.
```math
\begin{align*}
\overline{\text{total assets}} = \bar{B} + \bar{q} \bar{K}
\end{align*}
```

In steady state the absolute value of the total demanded assets of households $B^D_t$ (the sum of liquid asset holdings $< 0$) is equal to:

```math
\begin{align*}
\bar{B}^D = -\sum \left( \overline{\text{distr}}^b \cdot \mathbf{1}_{\{\text{grid}_b < 0\}} \cdot \text{grid}_b \right)
\end{align*}
```

#### Market clearing and further variables

This total value of demanded assets also enters the steady-state resource constraint:
```math
\begin{align*}
\bar{C} = \bar{Y} - \delta_0 \bar{K} - \bar{G} -  \bar{B}^D (\bar{r}^D - \bar{r}^L)
\end{align*}
```

The constants entering the household problem are fixed to their steady-state level. This applies to  $\tau^p$, $\tau^C$, $\sigma$, and $q$.

Moreover, the steady-state values of introduced lag variables are set equal to the steady-state values of the corresponding original variables. This applies to the lags of $B^{gov}_t$, $w^F_t$, $q_t$, and $\tau_t$.
