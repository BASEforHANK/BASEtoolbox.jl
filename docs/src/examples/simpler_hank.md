# Simpler HANK model

## Equations

### Shocks
```math
\begin{align}
\log \text{Rshock}_{t+1} &= \rho_{\text{Rshock}} \log \text{Rshock}_{t} \\
\log \text{A}_{t+1} &= \rho_{\text{A}} \log \text{A}_{t} \\
\log \sigma_{t+1} &= \rho_{s} \log \sigma_{t}
\end{align}
```

### Lags
```math
\begin{align}
\log B_{\text{gov,lag},t+1} &= \log B_{\text{gov},t} \\
\log I_{\text{lag},t+1} &= \log I_{t} \\
\log w_{\text{lag},t+1} &= \log w_{t} \\
\log \tau_{\text{avg,lag},t+1} &= \log \tau_{\text{avg},t} \\
\log \tau_{\text{prog,lag},t+1} &= \log \tau_{\text{prog},t} \\
\log q\Pi_{\text{lag},t+1} &= \log q\Pi_{t}
\end{align}
```

### Actual equations

Monetary policy rule
```math
\begin{align}
\log RB_{t+1} &= RB  + \theta_{\pi} \log \pi_t + \log \text{Rshock}_t
\end{align}
```

Deficit rule and government budget constraint
```math
\begin{align}
\log B_{\text{gov},t+1} &= \log B_{\text{gov},t} \\
\log G_{t} &= \log \left( B_{\text{gov},t+1} + T_t - \frac{RB_t}{\pi_t} B_{\text{gov},t} \right)
\end{align}
```

Price Phillips curve
```math
\begin{align}
    \log \pi_{t} - \pi = \kappa \left(mc_t - \frac{1}{\mu} \right) + \beta \left(( \log \pi_{t+1} - \pi) \frac{Y_{t+1}}{Y_t} \right)
\end{align}
```

Wage Phillips curve
```math
\begin{align}
    \log \pi_{w,t} - \pi_w = \kappa_w \left(mc_{w,t} - \frac{1}{\mu_w} \right) + \beta \left(( \log \pi_{w,t+1} - \pi_w) \frac{N_{t+1} w_{t+1}}{N_t w_t} \right)
\end{align}
```

Real wage inflation
```math
\begin{align}
    \log \frac{w_t}{w_{t-1}} = \log \frac{\pi_{w,t}}{\pi_t}
\end{align}
```

#### Firms: Intermediate goods producers

Intermediate goods producers minimize the costs of production $w_t N_t + (r_t - 1 + q_t \delta_0) K_t$ subject to the production function $Y_t = K_t^{\alpha} N_t^{1-\alpha}$. Hence, the first order conditions are

```math
\begin{align}
    r_t &= \alpha \text{mc}_t \left( \frac{K_t}{N_t} \right)^{\alpha-1} + 1 - q_t \delta_0 \\
    w_t &= (1-\alpha) \text{mc}_t \left( \frac{K_t}{N_t} \right)^{\alpha}
\end{align}
```

Union profits
```math
\begin{align}
    \log \Pi^U_t = \log \left( w_t N_t (1 - \text{mc}_{w,t} ) \right)
\end{align}
```

Firm profits
```math
\begin{align}
    \log \Pi^F_t = \log \left( Y_t (1-\text{mc}_t) + q_t \left(K_{t+1}- (1 - \delta_0) K_t \right) - I_t \right)
\end{align}
```

Value of profit shares
```math
\begin{align}
    \log \left( \frac{\text{RB}_{t+1}}{\pi_{t+1}} \right) = \log \left( \frac{(1-\iota) (q^\Pi_t - 1) + \omega \Pi^F_t}{q^\Pi_t - 1} \right)
\end{align}
```

Return on liquid assets
```math
\begin{align}
    \log \text{RL}_t = \log \left( \frac{\text{RB}_t B_{\text{gov},t} + \pi_t \left( (1-\iota) (q^\Pi_t - 1) + \omega \Pi^F_t \right)}{B_t} \right)
\end{align}
```

Total liquid assets
```math
\begin{align}
    \log B_t = \log \left( B_{\text{gov},t} + q^\Pi_t - 1 \right)
\end{align}
```

Capital price
```math
\begin{align}
    q_t = 1
\end{align}
```

Capital accumulation
```math
\begin{align}
    K_{t+1} = I_t + (1 - \delta_0) K_t
\end{align}
```
