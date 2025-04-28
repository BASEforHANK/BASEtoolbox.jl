
module Tools

using ..Types

using LinearAlgebra, CategoricalArrays, Roots, ForwardDiff, Distributions, Printf
using SpecialFunctions: erf
using FFTW: dct, ifft

export Brent,
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
    realpart,
    distrSummaries,
    gini

include("Tools/BrentsMethod.jl")
include("Tools/LinInterpols.jl")
include("Tools/LocateFcn.jl")
include("Tools/GCintegration.jl")
include("Tools/DCT.jl")
include("Tools/FastRoot.jl")
include("Tools/MarkovChain.jl")
include("Tools/SchurUtils.jl")
include("Tools/DualUtils.jl")
include("Tools/Pdf2cdf.jl")
include("Tools/Broyden.jl")
include("Tools/CentralDerivatives.jl")
include("Tools/DistrSummaries.jl")
end
