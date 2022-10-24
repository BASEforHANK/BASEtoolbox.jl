# Bayesian Solution and Estimation (BASE) for HANK toolbox.

Toolbox for estimating medium-scale HANK models based on "Shocks, Frictions, and Inequality in US Business Cycles" by Christian Bayer, Benjamin Born and Ralph Luetticke ([link to paper](https://www.benjaminborn.de/publication/bbl_inequality_2020/)).

Note that the toolbox is not a 1-for-1 replication package for the linked paper. In particular, the preset resolution is smaller and the estimation settings are tuned down to ensure a quick run, e.g. short MCMC chain and no computation of Hessian at mode (see the documentation for details).

There is also currently an issue with MKL on macOS. Therefore, we use OpenBLAS for macOS machines. Note, that MKL offers considerable speed-ups.

To get started, clone or download the repository. The full documentation can be found locally in `docs/build`.
