# Simpler HANK model with 2 sectors

In this example we stick with the simple HANK model but split up sectors into two distinct sectors, whereas nominal rigidities apply for only one of them. The model aggregates back to a Cobb-Douglas production function.

## Economics

The two sectors of the economy produce housing goods $H$ and service goods $S$, respectively.
While there is a continuum of service goods producers with differentiated service goods, the representative housing goods producer is takes the flexible housing prices as given.
Both products are further bundled into the final consumption good, which allows us to leave the household preferences over consumption and thus the household problem unchanged. We also stick with the government sector, asset market setup and implemented shocks of the [simpler HANK](#simpler_hank.md) example.

### Details

Housing goods are produced by a representative firm using capital:
```math
\begin{align*}
H_t = Z_t K_t \ . \tag{indexes.$H$}
\end{align*}
```
The service goods producers have market power over their differentiated services and face Calvo price rigidites. A service goods producer $f$ produces their service good $S(f)$ with only labor from labor packers as input:
```math
\begin{align*}
S_{t}(f) = Z_t N_{t}(f) . \tag{indexes.$S$}
\end{align*}
```
Both sectors produce with aggregate productivity $Z_t$.

There are labor unions with non-zero profits, who supply their labor to the competitive labor packer, as in the simpler HANK example.

Finally, there is a competitive representative firm who bundles the housing and service goods for final consumption with a Cobb-Douglas production function.

```math
\begin{align*}
    Y_t = H_t ^{\alpha} S_t ^{1-\alpha} \tag{indexes.$Y$}
\end{align*}
```

It hereby aggregates the differentiated labor services with a Dixit-Stiglitz aggregator:

```math
\begin{align*}
    S_t = \left(\smallint_0^1 S_t(f)^{\frac{\eta-1}{\eta}}df\right)^{\frac{\eta}{\eta-1}}
\end{align*}
```

Price rigidities now solely lie on the side of the service good producers. Thus the model's Phillips curve defines the gross inflation rate of the aggregate service price level $P^S_t$ instead of the final goods price.

```math
\begin{align*}
    \log \pi_{t}^S - \log \bar{\pi}^S = \kappa \left(mc_t - \frac{1}{\mu} \right) + \beta \left(( \log \pi_{t+1}^S - \log \bar{\pi}) \frac{S_{t+1}}{S_t} \right)
,
\tag{indexes.mc}
\end{align*}
```

This aggregate price level $P^S_t$ follows from cost minimization by the final goods producer, who demands the service goods, and is given by:

```math
\begin{align*}
    P^S_t = \left(\smallint_0^1 P^S_t(f)^{1-η}df\right)^{\frac{1}{1-η}}
\end{align*}
```
In the following we define $μ = η/(η − 1)$ as the target markup of the service goods producers.

Wages paid to labor unions are in real terms:
```math
\begin{align}
    w^F_t = \left(\frac{H_t}{S_t}\right)^{\alpha} mc_t^s (1-\alpha) Z_t \tag{indexes.$w^F$}
\end{align}
```

Plugging in the production functions of capital and labor for the housing and service goods shows that wages are effectively the same as in the simpler HANK version. Profits of the service goods providers however differ from the monopolistic firm profits inside a single sector and are given by:
```math
\begin{align*}
\Pi^S_t = (1-\alpha) Y_t(1- mc_t^s) \tag{indexes.$Π^S$}
\end{align*}
```

How these expressions for wages and profits can be derived is shown in the section below.

Final goods producers minimize the costs of production $P_t^S S_t + P^H_t H_t$ subject to the **Cobb Douglas** production function $Y_t = H_t^{\alpha} S_t^{1-\alpha}$. With zero profits, they charge the price

```math
\begin{align}
P^C_t = \left( \frac{ P^H_t}{\alpha} \right)^\alpha \left( \frac{P^S_t}{1-\alpha}\right)^{1-\alpha},
\end{align}
```
where $P^H$ and $P^S$ are the price of housing and the aggregate price level of service goods, respectively.

Following from equation (1) we can write the final consumption good's inflation as:
```math
\begin{align}
\pi_t = (\pi_{t}^S)^\alpha (\pi_{t}^K)^{1-\alpha} \tag{indexes.$\pi^S$}
\end{align}
```

The inflation rate of the price $P^H$ hereby is equal to the inflation of the rental factor of capital plus depreciation $\pi_{t}^K$, as the competitive housing producer cannot levy markups.


```math
\begin{align}
\pi_t^K = \frac{(RK_{t+1} + \delta_0)}{(RK_{t} + \delta_0)} \tag{indexes.$\pi^K$}
\end{align}
```

Profits and factor payments of the unions remain the same, as their provision of labor services to the labor packers doesn't change.

Formally, entrepreneur profits include both service and housing goods producer profits. However, housing firm profits are merely defined as an auxiliary variable, as they are zero in equilibrium.
```math
\begin{align*}
\Pi_{E,t} = \Pi_{S,t} + \Pi_{H,t} \tag{indexes.$Π_E$}
\end{align*}
```
```math
\begin{align*}
\text{with } \Pi_H = q_t \left( K_{t+1} - (1 - \delta_0) K_{t} \right) - I_t
\end{align*}
```
### Derivations of Equilibrium Conditions

!!! note
    Because there no longer is a single product, the aggregate production price of the final consumption good is denoted $P^c_t$ in the following derivations. In these derivations we drop the time index for simplicity.

#### Derivation of service goods producers paid wages (to unions)

Analogous to $mc_{t}^w$ of the unions in BBL (2024), $mc_{t}^S$ is the actual markdown of wages the service good providers pay to unions, so the following holds:

```math
\begin{align}
P^S &= \frac{1}{mc^s} \cdot \frac{W^F}{Z} \nonumber \\
P^S &= \frac{1}{mc^s} \cdot \frac{w^F}{Z} P^C  \quad \text{using real wages.}
\end{align}
```

They are also equal to the nominal marginal cost of the service providers devided by their obtained price, after rearranging:
```math
\begin{align*}
mc^s = \frac{1}{P^s} \cdot \frac{w^F P^C}{Z}
\end{align*}
```

Because housing demand is optimal, the price of housing is defined through its marginal contribution to the consumption bundle:

```math
\begin{align}
P^H = \alpha \frac{Y}{H}P^C
\end{align}
```

!!!tip
    The real price of housing $P^H/P^C$ in equilibrium will  then be equal to the rental rate of capital plus the cost of capital depreciation $(r_t + q_t δ_0)$.

Equivalently, optimal demand of service goods gives that:
```math
\begin{align}
P^S = (1-\alpha)\frac{Y}{S}P^C.
\end{align}
```

Plugging equations (2) and (3) for the prices of the two sector's goods into equation (1) for $P^C$, we get:

```math
\begin{align*}

P^C = \left( \frac{Y}{H} \right)^\alpha (P^C)^\alpha \left( \frac{1}{mc^s}\frac{w^F P^C }{Z (1-\alpha)} \right)^{1-\alpha}\\
\iff (P^C)^{1-\alpha} = \left( \frac{Y}{H} \right)^\alpha \left( \frac{1}{mc^s}\frac{w^F P^C }{Z (1-\alpha)} \right)^{1-\alpha}\\
P^C = \left( \frac{Y}{H} \right)^{\frac{\alpha}{1-\alpha}} \frac{1}{mc^s}\frac{w^F P^C }{Z (1-\alpha)}

\end{align*}
```

Finally this gives for real wages paid by firms that:


```math
\begin{align*}
w^F &= \left( \frac{H}{Y}\right)^{\frac{\alpha}{1-\alpha}} mc^s (1-\alpha) Z \nonumber \\
\iff  w^F &= \left(\frac{H}{S}\right)^\alpha mc^s (1-\alpha) Z
\end{align*}
```

Plugging in $P^S/P^C$ from equation (4) into equation (2) will yield the same result.

#### Derivation of service goods producers profits
Profits of the service goods producer are in real terms:

```math
\begin{align*}
\Pi^S = \frac{S (P^S-W^F/Z)}{P^C} = S(P^S/P^C - w^F/Z)
\end{align*}
```

Plugging in the relationship of wages and service goods prices of equation (2) following the definition of realized markdowns gives:

```math
\begin{align*}
\Pi^S = S ( P^S / P^C - P^S /P^C mc^s ) = S \cdot P^S/P^C  \cdot (1-mc^s)
\end{align*}
```

Recall equation (4), so that $P^S/P^C=(1-\alpha)Y/S$ and thus:
```math
\begin{align*}
\Pi^S = (1-\alpha) Y(1- mc^s)
\end{align*}
```

## Steady state

The steady state of this model economy is characterized by multiple conditions, of which many are similar to the conditions of the [simpler HANK](#simpler_hank.md) example.
All conditions that differ are related to the additional variables and changes in the production sector.

The steady state values of two additional variables for housing goods and service goods are fixed by their respective input levels in steady state and total factor productivity:
```math
\begin{align*}
\bar{H} &= \bar{K} \bar{Z} \\
\bar{S} &= \bar{N} \bar{Z}
\end{align*}
```

Output and wages are written as functions of housing goods and service goods in the steady state conditions of this model and therefore become:
```math
\begin{align*}
\bar{Y} &= bar{H}^\alpha \bar{S}^{1 - \alpha}\\
\bar{w}^F &= \bar{mc} (1 - \alpha) \, \bar{Z} \left( \frac{\bar{H}}{\bar{S}} \right)^{\alpha}
\end{align*}
```

Besides consumption price and wage inflation, the additional inflation measures of service goods $\bar{\pi}^S$ and of capital goods $\bar{\pi}^K$ are normalized to one in steady state.

Lastly, firm profits are lower than in the model with a single sector also in the steady state. They are only equal to the monopolistic profits of the service sector:
```math
\begin{align*}
\bar{\Pi}^S = (1-\alpha) \bar{Y}(1- \bar{mc}^s)
\end{align*}
```
