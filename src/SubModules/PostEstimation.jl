# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------

module PostEstimation

    # Sibling modules
    using ..PerturbationSolution

    # 3rd Party modules
    using   Plots,
            VegaLite,
            StatsPlots,
            LinearAlgebra,
            CategoricalArrays,
            DataFrames
            
    using FFTW: ifft

        
    include("PostEstimation/compute_hist_decomp.jl")
    include("PostEstimation/compute_irfs_vardecomp.jl")
    include("PostEstimation/compute_vardecomp_bounds.jl")
    include("PostEstimation/plot_irfs.jl")
    include("PostEstimation/plot_vardecomp.jl")
    include("PostEstimation/compute_bcfreq_vardecomp.jl")
    
    export compute_irfs_vardecomp,
        plot_irfs,
        compute_hist_decomp,
        plot_vardecomp,
        compute_bcfreq_vardecomp,
        compute_vardecomp_bounds

end # end submodule Estimation