# __precompile__(false)
# Code runs on Julia 1.9.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------
# Packages used: Plots Distributions BenchmarkTools JLD2 FileIO DataFrames ForwardDiff
# SparseArrays LinearAlgebra Random LaTeXStrings MatrixEquations Roots KrylovKit JSON 
# CodecZlib SpecialFunctions FFTW Parameters Setfield MCMCChains StatsPlots Optim CSV 
# OrderedCollections Flatten FieldMetadata MKL

module BASEforHANK

if !Sys.isapple() # issues encountered when using mkl with macos + more than 1 thread
    using MKL
end

using   LinearAlgebra,
        JLD2,
        FileIO,
        Setfield,
        Flatten,
        Distributions,
        Random,
        FieldMetadata

import  Flatten: flattenable
import  FieldMetadata: label
# Submodules only required by sibling modules
include("SubModules/Tools.jl")
include("SubModules/EconFunc.jl")

# Submodules that define functions used by parent
include("SubModules/Parsing.jl")
using .Parsing
include("SubModules/Steady.jl")
using .Steady
include("SubModules/Macro.jl")
using .Macro
include("SubModules/Estimation.jl")
using .Estimation
include("SubModules/PostEstimation.jl")
using .PostEstimation


# Structs to export
export  ModelParameters,
        NumericalParameters,
        EstimationSettings,
        SteadyResults,
        LinearResults,
        EstimResults,
        SteadyState

# Own Functions to export
export  compute_steadystate,
        call_find_steadystate,
        call_prepare_linearization,
        linearize_full_model,
        model_reduction,
        update_model,
        find_mode,
        montecarlo,
        compare_2_linearizations,
        reduction_quality,
        reduction_quality_seq,
        compute_irfs_vardecomp,
        plot_irfs,
        compute_hist_decomp,
        plot_vardecomp,
        compute_bcfreq_vardecomp,
        compute_vardecomp_bounds

# Functions passed through from 3rd Party packages
export  mode,
        metaflatten,
        prior,
        label,
        jldsave,
        @set!,
        @load

include("Model/input_aggregate_names.jl")
include("Preprocessor/prior.jl")
# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------



@doc raw"""
    call_findsteadystate()

Computes the steady state and fills the SteadyState struct -- without further steps of preparing the linearization.

# Returns
`struct` `SteadyResults`, containing returns of [`find_steadystate()`](@ref)
"""
function call_find_steadystate(m_par)
    #Calculate steady state capital stock
    KSS, VmSS, VkSS, distrSS, n_par, m_par = find_steadystate(m_par)

    return SteadyState(KSS, VmSS, VkSS, distrSS, n_par)
end

@doc raw"""
    call_prepare_linearization()

Runs the prepare linearization and fills the SteadyResults struct, sr.

# Returns
`struct` `SteadyResults`, containing returns of [`find_steadystate()`](@ref)
"""
function call_prepare_linearization(ss, m_par)

    # Prepare steadys state information for linearization
    XSS,
    XSSaggr,
    indexes,
    indexes_r,
    indexes_aggr,
    compressionIndexes,
    n_par,
    m_par,
    CDFSS,
    CDF_m,
    CDF_k,
    CDF_y,
    distrSS = prepare_linearization(ss.KSS, ss.VmSS, ss.VkSS, ss.distrSS, ss.n_par, m_par)

    println("Number of DCTs for Vm:")
    println(length(compressionIndexes[1]))

    println("Number of DCTs for Vk:")
    println(length(compressionIndexes[2]))

    println("Number of DCTs for COP:")
    println(length(compressionIndexes[3]))


    return SteadyResults(
        XSS,
        XSSaggr,
        indexes,
        indexes_r,
        indexes_aggr,
        compressionIndexes,
        n_par,
        m_par,
        CDFSS,
        CDF_m,
        CDF_k,
        CDF_y,
        distrSS,
        state_names,
        control_names,
    )
end

@doc raw"""
    compute_steadystate()

Compute steady state including the preparation for linearization

# Returns
`struct` `SteadyResults`, containing returns of [`find_steadystate()`](@ref)
"""
function compute_steadystate(m_par)
    #Calculate steady state capital stock
    ss = call_find_steadystate(m_par)

    sr = call_prepare_linearization(ss, m_par)

    return sr
end



@doc raw"""
    linearize_full_model()

Linearize the full model (i.e. including idiosyncratic states and controls) around the steady state, and solves
using [`LinearSolution()`](@ref).

# Returns
`struct` `LinearResults`, containing
- `A::Array{Float64,2}`,`B::Array{Float64,2}`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`]
    and `XPrime` [`A`]
- `State2Control::Array{Float64,2}`: observation equation
- `LOMstate::Array{Float64,2}`: state transition equation
"""
function linearize_full_model(sr::SteadyResults, m_par::ModelParameters)
    A = zeros(sr.n_par.ntotal, sr.n_par.ntotal)
    B = zeros(sr.n_par.ntotal, sr.n_par.ntotal)

    if sr.n_par.verbose
        println("Initial linearization")
    end
    State2Control, LOMstate, SolutionError, nk, A, B =
        LinearSolution(sr, m_par, A, B; estim = false)

    return LinearResults(State2Control, LOMstate, A, B, SolutionError, nk)
end


@doc raw"""
    find_mode(sr, lr)

Find parameter that maximizes likelihood of data given linearized model `lr`.

# Arguments
- `sr::SteadyResults`
- `lr::LinearResults`

# Returns
`struct` `EstimResults`, containing all returns of [`mode_finding()`](@ref)
"""
function find_mode(sr::SteadyResults, lr::LinearResults, m_par::ModelParameters)
    if sr.n_par.verbose
        println("Started mode finding. This might take a while...")
    end
    if e_set.mode_start_file == ""
        priors = collect(metaflatten(m_par, prior)) # model parameters
        if e_set.me_treatment != :fixed
            append!(priors, e_set.meas_error_distr)         # add the meas. error priors
        end
        par_start = mode.(priors)
    else
        @load e_set.mode_start_file par_final
        par_start = copy(par_final)
    end
    par_final,
    hessian_final,
    posterior_mode,
    meas_error,
    meas_error_std,
    parnames,
    Data,
    Data_missing,
    H_sel,
    priors,
    smoother_output,
    m_par,
    sr,
    lr = mode_finding(sr, lr, m_par, e_set, par_start)

    if sr.n_par.verbose
        println("Mode finding finished.")
    end

    lr = update_model(sr, lr, m_par)

    er = EstimResults(
        par_final,
        hessian_final,
        meas_error,
        meas_error_std,
        parnames,
        Data,
        Data_missing,
        H_sel,
        priors,
    )

    return er, posterior_mode, smoother_output, sr, lr, m_par
end



@doc raw"""
    montecarlo(mr,er;file=e_set.save_posterior_file)

Sample posterior of parameter vector with [`rwmh()`](@ref), take sample mean as
parameter estimate, and save all results in `file`.

# Arguments
- `sr::SteadyResults`
- `mr::LinearResults`
- `er::EstimResults`
"""
function montecarlo(
    sr::SteadyResults,
    lr::LinearResults,
    er::EstimResults,
    m_par::ModelParameters;
    file::String = e_set.save_posterior_file,
)
    hessian_sym = Symmetric(nearest_spd(inv(er.hessian_final)))
    if sr.n_par.verbose
        println("Started MCMC. This might take a while...")
    end
    if e_set.multi_chain_init == true
        init_draw, init_success =
            multi_chain_init(er.par_final, hessian_sym, sr, lr, er, m_par, e_set)

        par_final = init_draw
        if init_success == false
            error("Couldn't find initial value that produces posterior")
        end
    else
        par_final = copy(er.par_final)
    end

    draws_raw, posterior, accept_rate =
        rwmh(par_final, hessian_sym, sr, lr, er, m_par, e_set)

    ##
    parnames_ascii = collect(metaflatten(m_par, label))
    if e_set.me_treatment != :fixed
        for i in eachindex(e_set.meas_error_input)
            push!(parnames_ascii, string("sigma_me_", e_set.meas_error_input[i]))
        end
    end

    chn = Chains(
        reshape(
            draws_raw[e_set.burnin+1:end, :],
            (size(draws_raw[e_set.burnin+1:end, :])..., 1),
        ),
        [string(parnames_ascii[i]) for i = 1:length(parnames_ascii)],
    )
    chn_summary = summarize(chn)
    par_final = chn_summary[:, :mean]

    ##
    if e_set.me_treatment != :fixed
        m_par =
            Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(er.meas_error)])
    else
        m_par = Flatten.reconstruct(m_par, par_final)
    end

    lr = update_model(sr, lr, m_par)

    smoother_output = likeli(par_final, sr, lr, er, m_par, e_set; smoother = true)

    if sr.n_par.verbose
        println("MCMC finished.")
    end
    return sr,
    lr,
    er,
    m_par,
    draws_raw,
    posterior,
    accept_rate,
    par_final,
    hessian_sym,
    smoother_output
end

end # module BASEforHANK
