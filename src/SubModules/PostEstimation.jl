module PostEstimation

# Sibling modules
using ..PerturbationSolution
using ..Tools
using ..Parsing

# 3rd Party modules
using Plots,
    StatsPlots, LinearAlgebra, CategoricalArrays, DataFrames, AlgebraOfGraphics, Printf

using LaTeXStrings
using CairoMakie: BarPlot, Label
using FFTW: ifft
using KernelDensity: kde
using Interpolations

include("PostEstimation/compute_irfs.jl")
include("PostEstimation/compute_vardecomp.jl")
include("PostEstimation/compute_hist_decomp.jl")

include("PostEstimation/plot_irfs.jl")
include("PostEstimation/plot_irfs_cat.jl")
include("PostEstimation/plot_distributional_irfs.jl")
include("PostEstimation/plot_vardecomp.jl")
include("PostEstimation/plot_hist_decomp.jl")

export compute_irfs,
    compute_vardecomp,
    compute_vardecomp_bcfreq,
    compute_hist_decomp,
    plot_irfs,
    plot_irfs_cat,
    plot_vardecomp,
    plot_vardecomp_bcfreq,
    plot_distributional_irfs,
    plot_hist_decomp

end
