# Bayesian Solution and Estimation (BASE) for HANK toolbox.

Toolbox for estimating medium-scale HANK models based on "Shocks, Frictions, and Inequality in US Business Cycles" by Christian Bayer, Benjamin Born and Ralph Luetticke ([link to paper](docs/BayerBornLuetticke_2024_AER.pdf)).

To get started, clone or download the repository. The full documentation can be found [here](https://baseforhank.github.io/BASEtoolbox.jl/) or locally in `docs/build`. Note that the toolbox is not a 1-for-1 replication package for the linked paper. In particular, the preset resolution is smaller and the estimation settings are tuned down to ensure a quick run, e.g. short MCMC chain and no computation of Hessian at mode (see the [documentation](https://baseforhank.github.io/BASEtoolbox.jl/) for details).
