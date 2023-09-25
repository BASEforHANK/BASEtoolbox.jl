# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------
# Packages used: Plots Distributions BenchmarkTools JLD2 FileIO DataFrames ForwardDiff
# SparseArrays LinearAlgebra Random LaTeXStrings MatrixEquations Roots KrylovKit JSON 
# CodecZlib SpecialFunctions FFTW Parameters Setfield MCMCChains StatsPlots Optim CSV 
# OrderedCollections Flatten FieldMetadata MKL

module BfH_Tools

#3rd Party modules
using   LinearAlgebra,
        CategoricalArrays,
        Roots,
        ForwardDiff,
        Distributions
using   SpecialFunctions: erf
using   FFTW: dct, ifft

export  Brent,
        broyden,
        CustomBrent,
        centralderiv,
        mydctmx,
        uncompress,
        compress,
        select_ind,
        Fastroot,
        my_integrate,
        myinterpolate3,
        mylinearinterpolate,
        mylinearinterpolate2,
        mylinearinterpolate3,
        mylinearinterpolate_mult2,
        mylinearinterpolate_mult3,
        mylinearinterpolate!,
        mylinearinterpolate_mult2!,
        mylinearinterpolate_mult2!,
        mylinearinterpolate_mult3!,
        locate,
        Tauchen,
        ExTransition,
        cdf_to_pdf,
        pdf_to_cdf,
        real_schur,
        dualpart,
        realpart

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------

include("BfH_Tools/BrentsMethod.jl")
include("BfH_Tools/LinInterpols.jl")
include("BfH_Tools/LocateFcn.jl")
include("BfH_Tools/GCintegration.jl")
include("BfH_Tools/DCT.jl")
include("BfH_Tools/FastRoot.jl")
include("BfH_Tools/MarkovChain.jl")
include("BfH_Tools/Schur_and_DualUtils.jl")
include("BfH_Tools/Pdf2cdf.jl")
include("BfH_Tools/Broyden.jl")
include("BfH_Tools/centralderivatives.jl")

end # module BfH_Tools