# Post Estimation Commands

The `PostEstimation` submodule provides a set of functions for post-estimation analysis. These functions are designed to work with the output of the model estimation process and provide further insights into the model's behavior. The functions provide plots of impulse responses, variance decompositions and historical decompositions, as well as tables in form of dataframes, where appropriate.  

## Usage

You can use these functions by importing the `PostEstimation` module. Here's an example:

```julia
using .BASEforHANK.PostEstimation

# Now you can call the functions
compute_hist_decomp(...)
```

## Functions
### Impulse Responses
```@docs
compute_irfs_vardecomp
plot_irfs
```
### Variance Decomposition
```@docs
compute_bcfreq_vardecomp
plot_vardecomp
compute_vardecomp_bounds
```
### Historical Decomposition
```@docs
compute_hist_decomp
```


