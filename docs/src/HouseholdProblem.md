# Documentation of the household problem

<span style="color:red">**TO DO: Note: this is a preliminary version, needs to be checked carefully and written more clearly.**</span>

This section documents the household problem that is at the core of the package. Most features are also present in BBL, however, the package extends the model. As a user, given the household problem, you can specify the "surrounding" aggregate equations that determine the inputs to the household problem. You can switch off features of the household problem by setting the parameters or variables accordingly.

The **most important assumptions** of the household problem are the following:

1. GHH preferences with labor supply decision
2. Taxes: linear and progressive labor and entrepreneurial income tax, linear value-added tax, linear union-profit tax
3. Entrepreneurship as an idiocyncratic state
4. 2 accounts: liquid and illiquid assets

At the end of this section, we provide a list of **necessary inputs** to the household problem.

## Model setup

There is a continuum of ex-ante identical households of measure one, indexed by $i$. Households are infinitely lived, have time-separable preferences with time-discount factor $\beta$, and derive felicity from consumption $c_{it}$ and leisure. They obtain income from supplying labor, $n_{it}$, from renting out capital, $k_{it}$, and from earning interest on bonds, $b_{it}$. What is more, they (potentially) receive entrepreneurial profits, $\Pi_t^E$, and union profits, $\Pi_t^U$. Households pay taxes on labor and profit income as well as a value-added consumption tax.

The model features household heterogeneity. Households differ in their productivity and in whether they obtain profit income. They face incomplete markets, and capital as an asset is illiquid while bonds are liquid.

## Preferences

With respect to leisure and consumption, households have  GHH preferences and maximize the discounted sum of felicity
```math
\begin{align*}
    \mathbb{E}_0 \sum_{t=0}^\infty \beta^t u\left[c_{it} - G(h_{it},n_{it})\right].
\end{align*}
```
The maximization is subject to the budget constraints described further below. The felicity function $u$ exhibits a constant relative risk aversion (CRRA) with risk aversion parameter $\xi$,
```math
\begin{align*}
    u(x_{it}) &= \frac{x_{it}^{1-\xi}}{1-\xi}, \tag{CRRA}
\end{align*}
```
where $x_{it} = c_{it} - G(h_{it},n _{it})$ is household $i$'s composite demand for goods consumption, $c_{it}$, and leisure, and $G$ measures the disutility from work. While $n_{it}$ denotes household $i$'s labor supply, $h_{it}$ is the household's labor productivity.

$G(h_{it}, n_{it})$ measures the disutility from work. We assume the following form for the disutility from work with $\gamma>0$:
```math
\begin{align*}
    G(h_{it}, n_{it}) &= h_{it}^{1-\tau^p} \frac{n_{it}^{1+\gamma}}{1+\gamma}.
\end{align*}
```
This functional form has several advantages that we will use below. First, the constant elasticity of $G$ with respect to $n_{it}$ simplifies the expression for the composite consumption good, $x_{it}$, in the budget constraint. When the Frisch elasticity of labor supply $\gamma$ is constant and the tax schedule has the form specified below, the disutility of labor is always a fraction of labor income and constant across households. Therefore, in both the household's budget constraint and felicity function, only after-tax income enters and neither hours worked nor productivity appear separately. Second, scaling the disutility of labor by $h_{it}^{1−\tau^p}$ ensures a normalization to one in the stationary equilibrium (see BBL p. 1220/1221 for details). In the exponent thereof, $\tau^p$ represents the steady state level of progressivity of the tax code, explained below. The functional form simplifies the household problem as $h_{it}$ drops out from the first-order condition as long as tax progressivity is constant. Then, all households supply the same number of hours. Yet, we allow for tax progressivity to vary over time below. Still, with this functional form, the household problem is more tractable.

## Income taxation

We assume a progressive labor and entrepreneurial income tax schedule as in BBL (with one modification explained below) given by
```math
\begin{align*}
    f_\tau(x) = (1-\tau^l_t) x^{1-\tau^p_t}  (\bar{Y}^g_t)^{\tau^p_t}. \tag{Tax func}
\end{align*}
```
Here (and this is the modification relative to BBL), the average level of gross income across households, $\bar{Y}^g_t$, influences individual taxing and thus net income, $y^n_{it}$. $\tau^l_t$ and $\tau^p_t$ determine the level and progressivity of the tax code.

The tax level shifter $(\bar{Y}^g_t)^{\tau^p_t}$ ensures that in the corresponding representative agent economy, tax income (holding behavior constant) is independent of progressivity.
To see this, plug in the average gross income, $\bar{Y}^g_t$, which corresponds to the gross income of a representative agent, as tax basis $x$. This gives $f_\tau(\bar{Y}^g_t)= (1-\tau^l_t)\bar{Y}^g_t$, so that in the corresponding representative agent economy, the tax progressivity vanishes.
In the code, this tax level shifter can be easily switched off.

Income related to union profits, $\Pi^U_t$, is taxed at a rate $\bar \tau_t$ that is independent of the recipient's labor income.

Asset income is not taxed.

## Idiosyncratic productivity

The household sector is subdivided into two types of agents: workers and entrepreneurs. The transition between both types is stochastic. On top, workers face idiosyncratic labor productivity risk. Both, workers and entrepreneurs, rent out physical capital, but only workers supply labor. Entrepreneurs do not work, but earn entrepreneurial profits, $\Pi_t^E$. Profits of unions, $\Pi_t^U$, are equally distributed across workers. All households self-insure against the income risks they face by saving in a liquid asset (bonds) and a less liquid asset (capital). Trading these illiquid assets is subject to random participation in the capital market.

We assume that productivity evolves according to a log-AR(1) process and a fixed probability of transition between the worker and the entrepreneur state:
```math
\tilde{h}_{it} =
\begin{cases}
\exp\left(\rho_{h} \log \tilde{h}_{it-1} + \epsilon^{h}_{it}\right) & \text{with probability } 1-\zeta \text{ if } h_{it-1} \neq 0, \\
1 & \text{with probability } \iota \text{ if } h_{it-1} = 0, \\
0 & \text{else}.
\end{cases}
```
with individual productivity $h_{it}= \frac{\tilde{h}_{it}}{\tilde{H}_{it}}$ such that $\tilde h_{it}$ is scaled by its cross-sectional average, $\tilde{H}_{it} = \int \tilde h_{it} di$, to make  sure that average worker productivity is constant. The shocks $\epsilon^{h}_{it}$ to productivity $\tilde{h}_{it}$ are normally distributed with variance $\bar \sigma^2_{h,t}$.

With probability $\zeta$ households become entrepreneurs ($h=0$). With probability $\iota$ an entrepreneur returns to the labor force with median productivity. Besides their labor income, workers receive a share in union rents, $\Pi^U_t$, which are distributed lump sum, leading to labor-income compression. For tractability, we assume union profits to be taxed at a rate $\bar{τ}_t$ independent of the recipient's labor income.

## Budget constraint

A household's net labor income, $y^n_{it}$, is given by
```math
\begin{align*}
    y^n_{it} = f_\tau (w^H_t h_{it} n_{it})
\end{align*}
```
where $w^H_t$ is the real wage rate that household receive, $h_{it}$ is the (scaled) productivity of the household, and $n_{it}$ is the amount of labor supplied.

Given incomes, households optimize intertemporally subject to their budget constraint
```math
\begin{align*}
(1 + \tau^c_t) c_{it} + b_{it+1} + q_t k_{it+1} &= RR^{i}_{t} b_{it} + (q_t + r^K_t)k_{it} + y^n_{it} + 1\{h_{it} \neq 0\}(1 - \bar \tau_t) \Pi^U_t + 1\{h_{it} = 0\} f_\tau(\Pi^E_t),
\end{align*}
```
where $\Pi^U_t$ is union profits taxed at the average tax rate of labor income $\bar \tau_t$, $\Pi^E_t$ is profit payouts to entrepreneurs, $b_{it}$ is real liquid assets, $k_{it}$ is the amount of illiquid assets, $q_{t}$ is the price of these assets, $r^K_{t}$ is their dividend, and $RR^{i}_{t}$ is the gross real interest rate on liquid assets, which depends on whether the household (thus, the dependence on $i$) is a borrower or lender.

All households that do not participate in the capital market $(k_{it+1}=k_{it})$ still obtain dividends and can adjust their liquid asset holdings. Depreciated capital has to be replaced for maintenance, such that the dividend, $r^K_{t}$, is the net return on capital. Holdings of bonds have to be above an exogenous debt limit $\underline{B},$ and holdings of capital have to be non-negative, that is
```math
k_{it+1} \geq 0, \quad b_{it+1} \geq \underline{B}.
```

Recall, $f_\tau(\Pi_t^E)$ is the progressive tax function defined above, applied to profit payouts to entrepreneurs. In contrast, merely the average tax rate $\bar \tau_t$ (further details below) is applied to the union's profits.

Additionally, consumption is taxed at a rate $\tau^c_t$, i.e. a VAT tax.

### Side note: computation of the budget constraint elements in the code

In the code, the elements of the budget constraint that are necessary to compute the household problem are computed by the `incomes()` function. Based on the parameters and aggregate variables, the function computes (at least) the following elements:
1. Net labor and union income for workers, adjusted for the composite good / net entrepreneurial profits for entrepreneurs, that is $y^n_{it} + 1\{h_{it} \neq 0\}(1 - \bar \tau_t) \Pi^U_t + 1\{h_{it} = 0\} f_\tau(\Pi^E_t)$ plus the adjustment for the composite good derived below.
2. Rental income from illiquid assets, that is $r^K_t k_{it}$.
3. Liquid asset income, that is $RR^{i}_{t} b_{it}$.
4. Liquidation income from illiquid assets, that is $q_t k_{it+1}$.

## Asset market

Households make their savings and portfolio choice between liquid bonds and illiquid capital in light of a capital market friction that renders capital illiquid because participation in the capital market is random and i.i.d. in the sense that only a fraction $\lambda$ of households are selected to be able to adjust their capital holdings in a given period.

The ex-post real return on the liquid asset, $RR^{i}_{t}$, is
```math
\begin{align*}
    RR^{i}_{t} =
\begin{cases}
    RR^{L}_{t} & \text{ if } b_{it} \geq 0, \\
    RR^{D}_{t} & \text{ else}. \tag{Return liquid}
\end{cases}
\end{align*}
```
Here, $RR^{L}_{t}$ is the real lending rate and $RR^{D}_{t}$ is the real borrowing rate. There is potentially a wedge between the two rates.

## Underlying decision problem

Households face the following decision problem:
```math
\begin{align*}
    \max_{\left\{c_{it},n_{it},b_{it+1},k_{it+1}\right\}} \mathbb{E}_0 \sum_{t=0}^\infty \beta^t u\left[c_{it} - G(h_{it},n_{it})\right],
\end{align*}
```
subject to the budget constraint and both short-selling constraints.

Moreover, as stated above, the household might not be able to adjust its capital holdings in a given period.
Expectations are thus formed over income states and the stochastic participation on the capital market.
In periods, where the household cannot adjust, it faces the same problem except for the additional condition $k_{it+1} = k_{it}$.

## Labor supply decision

Consider the following Lagrangian of the household problem (ignoring the asset choice here, since it is not relevant for the labor supply decision):
```math
\begin{align*}
\mathcal{L} &= E_0 \sum_{t=0}^\infty \beta^t \Bigg( u(c_{it} - G(h_{it}, n_{it})) + \lambda_{it} \left[ RR^{i}_{t} b_{it} + (q_t + r^K_t) k_{it} + y^n_{it} + 1\{h_{it} \neq 0\}(1 - \bar \tau_t)\Pi^U_t + 1\{h_{it} = 0\} f_\tau(\Pi^E_t) - (1 + \tau^c_t) c_{it}
 - b_{it+1} - q_t k_{it+1}\right] + \dots \Bigg)
\end{align*}
```
Remember, $y^n_{it}$ denotes the household's net labor income.

Combining the first-order conditions with respect to labor $n_{it}$ and consumption $c_{it}$ yields:
```math
\begin{align*}
\frac{\partial G(h_{it}, n_{it})}{\partial n_{it}} &= \frac{1}{(1 + \tau^c_t)}\frac{\partial y^n_{it}}{\partial n_{it}} = \frac{1}{(1 + \tau^c_t)} (1- \tau^p) \frac{y^n_{it}}{n_{it}}. \tag{FOC labor}
\end{align*}
```
Together with $\partial G(h_{it}, n_{it})/ \partial n_{it} = (1+\gamma)[G(h_{it}, n_{it})/n_{it}]$ (see above), we can express the composite consumption good, $x_{it}$, as
```math
\begin{align*}
x_{it} = c_{it} - G(h_{it}, n_{it}) &= c_{it} - \frac{1}{1+\gamma} \frac{1-\tau^p_t}{1+\tau^c_t}  y^n_{it}. \tag{Cons comp labor}
\end{align*}
```

We allow for the parameter $\tau^p_t$ governing the progressivity of the tax schedule to vary over time. In this case, as explained in footnote 12 of BBL, whenever tax progressivity does not coincide with its steady state level $\tau^p$, individual hours worked differ across agents. The labor supply under time-varying tax progressivity can be derived from eq. (FOC labor) and is similar to BBL equation (19a). It is given by:
```math
n_{it} = \left[ \frac{(1 - \tau^p_t)(1 - \tau^l_t)}{(1+\tau^c_t)}  (\bar{Y}^g_t)^{\tau^p_t} \right]^{\frac{1}{\gamma + \tau^p_t}}
(w^H_t)^{\frac{1 - \tau^p_t}{\gamma + \tau^p_t}}
h_{it}^{\frac{\tau^p - \tau^p_t}{\gamma + \tau^p_t}}
```

On aggregate, this results in total effective hours (compare equation (19b) of BBL):
```math
N_t = \int n_{it}h_{it} \, di = \left[ \frac{(1 - \tau^p_t)(1 - \tau^l_t)}{(1+\tau^c_t)}  (\bar{Y}^g_t)^{\tau^p_t} \right]^{\frac{1}{\gamma + \tau^p_t}}
(w^H_t)^{\frac{1 - \tau^p_t}{\gamma + \tau^p_t}} \underbrace{\int h_{it}^{\frac{\gamma + \tau^p}{\gamma + \tau^p_t}} \, di}_{=: H^p_t}, \tag{Total hours}
```

Individual hours worked, $n_{it}$, can thus be written as a function of total effective hours $N_{t}$: $n_{it}  = N_t/H^p_t  h_{it}^{(\tau^p - \tau^p_t)/(\gamma + \tau^p_t)}$.

Given the labor decision, net income can be calculated.

We plug the optimal supply of hours, $n_{it}$, into the equation for net income, $y^n_{it}$, which is then (compared to equation (22) in BBL):
```math
\begin{align*}
y^n_{it} &= (1 - \tau^l_t)(w^H_t h_{it} n_{it})^{1 - \tau^p_t} (\bar{Y}^g_t)^{\tau^p_t} = (1 - \tau^l_t)^{\frac{1+\gamma}{\gamma + \tau^p_t}}
 (1 - \tau^p_t) ^ {\frac{1-\tau^p_t}{\gamma + \tau^p_t}}
 (1 + \tau^c_t) ^ {-\frac{1-\tau^p_t}{\gamma + \tau^p_t}}
 w_t ^ {\frac{(1 + \gamma)(1 - \tau^p_t)}{\gamma + \tau^p_t}}
 h_{it} ^ {\frac{(\gamma + \tau^p)(1-\tau^p_t)}{\gamma + \tau^p_t}}
 (\bar{Y}^g_t)^{\frac{\tau^p_t (1+\gamma)}{\gamma + \tau^p_t}}.
 \end{align*}
```

Plugging individual hours as a function of total effective hours into gross income, results in
```math
y^g_{it} = w^H_t h_{it} n_{it} = w^H_t  \frac{N_t}{H^p_t}  h_{it}^\frac{\gamma + \tau^p}{\gamma + \tau^p_t}. \tag{Gross income}
```

Finally, consumption in the budget constraint can be replaced with the composite consumption good $x_{it}$:
```math
\begin{align*}
(1 + \tau^c_t) x_{it} + b_{it+1} + q_t k_{it+1} &= RR^{i}_{t} b_{it} + (q_t + r^K_t)k_{it} + \left(1 - \frac{1-\tau^p_t}{1+\gamma}\right) y^n_{it} + 1\{h_{it} \neq 0\}(1 - \bar \tau_t) \Pi^U_t + 1\{h_{it} = 0\} f_\tau(\Pi^E_t), \tag{BC with x}
\end{align*}
```
where $\left(1-\frac{1-\tau^p_t}{1+\gamma}\right)$ simplifies to $\left(\frac{\tau^p_t + \gamma}{1+\gamma}\right)$. This GHH factor is used in the code to transform net labor income, $y^n_{it}$, into the composite consumption good units. Here, we have to take the VAT additionally into account.

## Portfolio choice

### Description of the dynamic problem

In the following, unlike in other sections of the documentation, we will adopt recursive notation (where a prime, `'`, denotes next period's value) for clarity. Policy functions are denoted with a star, `*`, and subscripts indicate whether the policy pertains to the adjustment case, `a`, or the non-adjustment case, `n`. For example, $x^*_a$ represents the policy function for the consumption-leisure bundle in the adjustment case, while $x^*_n$ represents it in the non-adjustment case. For notational brevity, we refer to the consumption-leisure composite simply as consumption throughout.

Since a household's saving decision will be some non-linear function of that household's wealth and productivity, inflation and all other prices will be functions of the joint distribution, $\Theta_t$, of $(b,k,h)$ in $t$. This makes $\Theta_t$ a state variable of the household's planning problem. It evolves as a result of the economy's reaction to aggregate shocks according to:

```math
\begin{align*}
    \Theta_{t+1}(b',k',h') &=\lambda \int _{b'=b_{a,t}^*(b,k,h),k'=k_t^*(b,k,h)}\Pi(h,h')d\Theta_{t}(b,k,h)
    + (1-\lambda) \int_{b'=b_{n,t}^*(b,k,h),k'=k}\Pi(h,h')d\Theta_{t}(b,k,h)\ ,
\end{align*}
```
where $\Pi(\cdot)$ is the transition probability for $h$ and $b_{a/n,t}^*$ and $k_t^*$ are the time-$t$ optimal policies.

For simplicity, we summarize all effects of aggregate state variables, including the distribution of wealth and income, by writing the dynamic planning problem with time-dependent continuation values.

This leaves us with two Bellman equations, that characterize the household's problem: $V^{a}_t$ for the case where the household adjusts its capital holdings and $V^{n}_t$ for the case in which it does not adjust. A third function defines the expected continuation value, $\mathbb{W}_t$, over both,
```math
\begin{align*}
V^{a}_t(b,k,h) =& \max_{b_{a}',k'}  u[x(b,b_{a}',k,k',h)]+ \beta\mathbb{E}_t \mathbb{W}_{t+1}(b_{a}',k',h')\ , \nonumber\\
V^{n}_t(b,k,h) =& \max_{b_{n}'}  u[x(b,b_{n}',k,k,h)] + \beta \mathbb{E}_t\mathbb{W}_{t+1}(b_{n}',k,h')\ , \\
\mathbb{W}_{t+1}(b',k',h') =& \lambda V_{t+1}^{a}(b',k',h') +(1-\lambda)V^{n}_{t+1}(b',k',h') \ .\nonumber
\end{align*}
```
Expectations about the continuation value are taken with respect to all stochastic processes (i.e. with respect to idiosyncratic and aggregate variables) conditional on the current states.

The budget sets for the adjustment case $\Gamma_a$ and the non-adjustment case $\Gamma_n$ are given by
```math
\begin{align*}
    \Gamma_{a}(b,k,h) &= \Bigg\{ k' \geq 0, b' \geq \underline{B} \mid q (k' - k) + b' \leq \left(\frac{\tau^p + \gamma}{1+\gamma}\right) y^n + RR(b) b + r^K k \\
    & + 1\{h \neq 0\}(1- \tau)\Pi^U + 1\{h = 0\} f_\tau(\Pi^E) \Bigg\} \\
    \Gamma_{n}(b,k,h) &= \Bigg\{b' \geq \underline{B} \mid  b'\leq \left(\frac{\tau^p + \gamma}{1+\gamma}\right) y^n + RR(b) b + r^K k \\
    & + 1\{h_{it} \neq 0\}(1- \tau)\Pi^U + 1\{h_{it} = 0\} f_\tau(\Pi^E) \Bigg\}.
\end{align*}
```
Moreover, the consumption leisure composite $x$ is
```math
\begin{align*}
    x(b,b',k,k',h) &= \frac{1}{1+\tau^c} \Bigg[ \left(\frac{\tau^p + \gamma}{1+\gamma}\right) y^n + RR(b) b + r^K k \\
    & + 1\{h \neq 0\}(1- \tau)\Pi^U + 1\{h = 0\} f_\tau(\Pi^E) - b' - q (k' - k) \Bigg].
\end{align*}
```

As [Bayer, Luetticke, Pham-Dao and Tjaden (2019)](https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA13601) show in the online Appendix B, the Bellman equation fullfills Blackwell's sufficient condition for a contraction on the set of bounded, continuous and weakly concave functions, such that we can expect to find a solution to the problem.

When solving the problem, we do not use the value functions. Instead, we rely on the first-order conditions that characterize the optimal household behavior to solve for policy functions. The next section develops the equations that are necessary for an iterative algorithm that solves for optimal policies, before we describe the exact algorithm thereafter.

### Euler Equations


Denote the optimal policies for consumption, for liquid asset holdings and illiquid asset holdings as $x_a^{*}$, $b_a^{*}$, $k_a^{*}$, $x_n^{*}$, and $b_n^{*}$, respectively. The first-order conditions (FOCs) for a solution in the adjustment case and the non-adjustment case read:
```math
\begin{align*}
    k^{*}_a: &\quad  \frac{\partial u(x^{*}_{a})}{\partial x} \geq (1 + \tau^c) q^{-1} \beta \, \mathbb E \, \frac{\partial \mathbb W (b^{*}_{a} ,k^{*}_a ,h')}{\partial k} \tag{FOC1} \\
    b^{*}_{a}: &\quad \frac{\partial u(x^{*}_{a})}{\partial x} \geq (1 + \tau^c) \beta \, \mathbb E \, \frac{\partial \mathbb W (b^{*}_{a}, k^{*}_a, h')}{\partial b} \tag{FOC2} \\
    b^{*}_{n}: &\quad \frac{\partial u(x^{*}_{n})}{\partial x} \geq (1 + \tau^c) \beta \, \mathbb E \, \frac{\partial \mathbb W (b^{*}_{n}, k, h')}{\partial b} \tag{FOC3}
\end{align*}
```
Note the subtle difference between (FOC2) and (FOC3), which lies in the different illiquid assets stocks $k^{*}_a$ vs. $k$ entering as arguments on the right-hand-side. The difference arises from the asset market friction that only allows a fraction of households to adjust their illiquid asset holding. The continuation values (CVs) are defined as:
```math
\begin{align*}
    \frac{\partial \mathbb W}{\partial k} (b', k', h') &=  \lambda \frac{\partial V_{a}(b', k', h')}{\partial k} + (1 - \lambda)\frac{\partial V_{n}(b', k', h')}{\partial k} \tag{CV1} \\
    \frac{\partial \mathbb W}{\partial b} (b', k', h') &= \lambda \frac{\partial V_{a}(b', k', h')}{\partial b} + (1 - \lambda)\frac{\partial V_{n}(b', k', h')}{\partial b} \tag{CV2}
\end{align*}
```
To characterize continuation values, we calculate the derivatives of the value functions for the adjustment and non-adjustment cases with respect to liquid and illiquid assets. Using the Envelope conditions (ECs), we obtain:
```math
\begin{align*}
    \frac{\partial V_{a}}{\partial k}(b, k, h) &= \frac{\partial u[x^{*}_{a}(b, k, h)]}{\partial x} \left( \frac{q + r^K}{1 + \tau^c} \right) \tag{EC1} \\
    \frac{\partial V_{a}}{\partial b}(b, k, h) &= \frac{\partial u[x^{*}_{a}(b, k, h)]}{\partial x} RR(b) \frac{1}{1 + \tau^c} \tag{EC2} \\
    \frac{\partial V_{n}}{\partial b}(b, k, h) &= \frac{\partial u[x^{*}_{n}(b, k, h)]}{\partial x} RR(b) \frac{1}{1 + \tau^c} \tag{EC3} \\
    \frac{\partial V_{n}}{\partial k}(b, k, h) &= \frac{\partial u[x^{*}_{n}(b, k, h)]}{\partial x} \left(\frac{r^K}{1 + \tau^c}\right) \\
    & \quad + \beta \, \mathbb E \left[ \lambda \frac{\partial V_{a}[b^{*}_{n}(b, k, h), k, h']}{\partial k} + (1 - \lambda)\frac{\partial V_{n}[b^{*}_{n}(b, k, h), k, h']}{\partial k}\right] \\
    &= \frac{\partial u[x^{*}_{n}(b, k, h)]}{\partial x} \left(\frac{r^K}{1 + \tau^c}\right) + \beta \, \mathbb{E} \frac{\partial \mathbb W [b^{*}_{n}(b, k, h), k, h']}{\partial k} \tag{EC4}
\end{align*}
```
such that the marginal value of the illiquid asset in the non-adjustment case is defined recursively.

By substituting the Envelope conditions (ECs) into the continuation values (CVs) and the composite expression into the first-order conditions (FOCs), we derive the Euler equations (EEs) that characterize optimal household behavior. Using a slightly simplified notation, these are given by the following equations:
```math
\begin{align*}
    \frac{\partial u[x_{a}^{*}(b, k, h)]}{\partial x} &\geq q^{-1} \beta \, \mathbb E (1 + \tau^c) \left[ \lambda  \frac{\partial u[x^{*}_{a}(b^{*}_{a}, k^{*}_a, h')]}{\partial x} \left( \frac{q'+{r^K}'}{1 + {\tau^c}'} \right) + (1-\lambda)\frac{\partial V_{n}}{\partial k} (b_{a}^{*},k^{*}_a, h') \right] \tag{EE1} \\
    \frac{\partial u[x_{a}^{*}(b,k, h)]}{\partial x} &\geq \beta \, \mathbb E RR(b^{*})' \left( \frac{1 + \tau^c}{1 + {\tau^c}'} \right) \left[ \lambda  \frac{\partial u[x^{*}_{a}(b^{*}_{a},k^{*}_a, h')]}{\partial x} + (1-\lambda)\frac{\partial u[x^{*}_{n}(b_{a}^{*},k^{*}_a, h')]}{\partial x}\right] \tag{EE2} \\
    \frac{\partial u[x_{n}^{*}(b,k, h)]}{\partial x}  &\geq \beta \, \mathbb E RR(b^{*})' \left( \frac{1 + \tau^c}{1 + {\tau^c}'} \right) \left[ \lambda  \frac{\partial u[x_{a}^{*}(b_{n}^{*},k, h')]}{\partial x} + (1-\lambda)\frac{\partial u[x_{n}^{*}(b_{n}^{*},k, h')]}{\partial x}\right] \tag{EE3}
\end{align*}
```
In the non-adjustment case, the intuition behind the Euler equation is standard: the household equates the marginal loss of utility today from saving an additional unit with the expected marginal utility gain from saving that unit.

In the adjustment case, the household also faces a portfolio problem. Intuitively, it selects the optimal combination of the liquid and illiquid assets by comparing the expected one-period return difference between the two assets, $\mathbb E \left( RR(b^{*})' - \frac{{r^K}' + q'}{q} \right)$, weighted by the probability-adjusted marginal utilities. (We abstract from the non-differentiability at $b=0$.) After having decided the optimal portfolio combination, the household then faces the identical tradeoff between consumption today and consumption tomorrow as in the non-adjustment case.

### Solving the dynamic planning problem

The endogenous grid algorithm (EGM) is used to solve the dynamic household problem with two assets.
The first-order conditions and Euler equations, which were derived above, form the foundation of the iterative algorithm.
The document [`Computational Notes.md`](Computational Notes.md) provides a step-by-step explanation of the computational procedure implemented in the code.

## Inputs to the household problem

### Inputs: Parameters

- time-discount factor, $\beta$
- risk aversion parameter, $\xi$
- Frisch elasticity of labor supply, $\gamma$
- steady-state progressivity of the tax code, $\tau^p$
- persistence of the productivity shock, $\rho_h$
- variance of the productivity shock, $\bar \sigma^2_{h}$
- probability of becoming an entrepreneur, $\zeta$
- probability of returning to the labor force, $\iota$
- debt limit, $\underline{B}$
- probability of adjusting illiquid assets, $\lambda$

### Inputs: Aggregate variables

- real wage rate that households receive, $w^H_t$
- level of the tax code, see tax function (Tax func), $\tau^l_t$
- progressivity of the tax code, see tax function (Tax func), $\tau^p_t$
- average gross income across households, see tax function (Tax func), $\bar{Y}^g_t$
- variance of the productivity shock, $\bar \sigma^2_{h,t}$
- union profits, $\Pi^U_t$
- entrepreneurial profits, $\Pi^E_t$
- tax rate on union profits, $\bar \tau_t$
- VAT tax rate, $\tau^c_t$
- price of illiquid assets, $q_t$
- real dividend on illiquid assets, $r^K_t$
- real lending rate, $RR^{L}_t$
- real borrowing rate, $RR^{D}_t$
- cross-sectional average of productivity before scaling, $\tilde{H}_t$
- cross-sectional average of productivity, weighted by tax progressivity factor, $H^p_t$
- total effective hours, $N_t$
