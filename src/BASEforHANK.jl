# __precompile__(false)

module BASEforHANK
if !Sys.isapple() # issues encountered when using mkl with macos + more than 1 thread
    using MKL
end

using LinearAlgebra,
    JLD2,
    FileIO,
    JSON,
    Setfield,
    Flatten,
    Distributions,
    Random,
    FieldMetadata,
    MCMCChains,
    Printf

# Documentation: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    @printf "Paths are not defined, assuming documentation mode. If you want to run the code, please define the paths to the model.\n"
end

import Flatten: flattenable

# Submodules only required by sibling modules
include("SubModules/Types.jl")
using .Types
include("SubModules/Tools.jl")

# Submodules that define functions used by parent
include("SubModules/Parsing.jl")
using .Parsing
include("SubModules/IncomesETC.jl")
include("SubModules/SteadyState.jl")
using .SteadyState
include("SubModules/PerturbationSolution.jl")
using .PerturbationSolution
include("SubModules/Estimation.jl")
using .Estimation
include("SubModules/PostEstimation.jl")
using .PostEstimation

# Structs to export
export ModelParameters,
    NumericalParameters,
    EstimationSettings,
    SteadyResults,
    LinearResults,
    EstimResults,
    SteadyState,
    AbstractMacroModel,
    OneAsset,
    TwoAsset,
    CompleteMarkets

# Own Functions to export
export compute_steadystate,
    call_find_steadystate,
    call_prepare_linearization,
    linearize_full_model,
    model_reduction,
    update_model,
    find_mode,
    sample_posterior,
    find_field_with_value

# Own functions to export: computations and plotting
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

# Functions passed through from 3rd Party packages
export mode, metaflatten, prior, label, jldsave, @set!, @load

# Documentation: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    include("../examples/baseline/Model/input_aggregate_names.jl")
else
    include(Main.paths["src_example"] * "/Model/input_aggregate_names.jl")
end
include("Preprocessor/prior.jl")

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------

@doc raw"""
    call_find_steadystate(m_par)

Computes the steady state of the household problem and fills the SteadyStateStruct struct
(without further steps of preparing the linearization).

# Arguments
- `m_par::ModelParameters`

# Returns
- `SteadyStateStruct`, containing returns of [`find_steadystate()`](@ref)
"""
function call_find_steadystate(m_par::ModelParameters)
    @printf "\n"
    @printf "Compute the steady state...\n"

    KSS, WbSS, WkSS, distrSS, n_par, m_par = find_steadystate(m_par)

    @printf "Compute the steady state... Done.\n"

    return SteadyStateStruct(KSS, WbSS, WkSS, distrSS, n_par)
end

@doc raw"""
    call_prepare_linearization(ss, m_par)

Prepares linearization and fills the SteadyResults struct.

# Arguments
- `ss::SteadyStateStruct`: Output of [`call_find_steadystate()`](@ref)
- `m_par::ModelParameters`

# Returns
- `SteadyResults`, containing returns of [`prepare_linearization()`](@ref)
"""
function call_prepare_linearization(ss::SteadyStateStruct, m_par::ModelParameters)
    @printf "\n"
    @printf "Prepare the linearization...\n"

    # Prepare steady state information for linearization
    XSS,
    XSSaggr,
    indexes,
    indexes_aggr,
    compressionIndexes,
    n_par,
    m_par,
    CDFSS,
    CDF_bSS,
    CDF_kSS,
    CDF_hSS,
    distrSS = prepare_linearization(ss.KSS, ss.WbSS, ss.WkSS, ss.distrSS, ss.n_par, m_par)

    # Check that the steady state is well-defined
    idx_not = findall(isnan.(XSS) .| isinf.(XSS))
    if !isempty(idx_not)
        @warn "Steady state is not well-defined: NaN or Inf in XSS!"
        for i in idx_not
            myfield = find_field_with_value(indexes, i, true)
            @printf "Variable: %s Value: %f\n" myfield XSS[i]
        end
    end

    if n_par.verbose
        @printf "Number of DCTs for Wb: %d\n" length(compressionIndexes[1])
        @printf "Number of DCTs for Wk: %d\n" length(compressionIndexes[2])
        @printf "Number of DCTs for COP: %d\n" length(compressionIndexes[3])
    end

    @printf "Prepare the linearization... Done.\n"

    return SteadyResults(
        XSS,
        XSSaggr,
        indexes,
        indexes,
        indexes_aggr,
        compressionIndexes,
        n_par,
        m_par,
        CDFSS,
        CDF_bSS,
        CDF_kSS,
        CDF_hSS,
        distrSS,
        state_names,
        control_names,
    )
end

@doc raw"""
    compute_steadystate(m_par)

A wrapper for [`call_find_steadystate()`](@ref) and [`call_prepare_linearization()`](@ref).

# Arguments
- `m_par::ModelParameters`

# Returns
- `SteadyResults`, containing returns of [`prepare_linearization()`](@ref)
"""
function compute_steadystate(m_par::ModelParameters)
    ss = call_find_steadystate(m_par)
    sr = call_prepare_linearization(ss, m_par)
    return sr
end

@doc raw"""
    linearize_full_model(sr, m_par; ss_only = false)

Linearize the full model (i.e. including idiosyncratic states and controls) around the
steady state, and solves using [`LinearSolution()`](@ref).

# Arguments
- `sr::SteadyResults`: Output of [`call_prepare_linearization()`](@ref)
- `m_par::ModelParameters`
- `ss_only::Bool`: If `true`, only the steady state is checked

# Returns
- `LinearResults`, containing

    - `State2Control::Array{Float64,2}`: Matrix of observation equation
    - `LOMstate::Array{Float64,2}`: Matrix of state transition equation
    - `A::Array{Float64,2}`,`B::Array{Float64,2}`: First derivatives of
      [`PerturbationSolution.Fsys()`](@ref) with respect to arguments `X` [`B`] and `XPrime`
      [`A`]
"""
function linearize_full_model(sr::SteadyResults, m_par::ModelParameters; ss_only = false)
    @printf "\n"
    @printf "Linearizing the full model...\n"

    A = zeros(sr.n_par.ntotal, sr.n_par.ntotal)
    B = zeros(sr.n_par.ntotal, sr.n_par.ntotal)

    if ss_only
        @printf "Linearizing the full model... SS only.\n"
        return LinearSolution(sr, m_par, A, B; ss_only = true)
    end

    State2Control, LOMstate, SolutionError, nk, A, B =
        LinearSolution(sr, m_par, A, B; allow_approx_sol = true)

    @printf "Linearizing the full model... Done.\n"

    return LinearResults(State2Control, LOMstate, A, B, SolutionError, nk)
end

@doc raw"""
    find_mode(sr, lr, m_par, e_set)

Find parameter that maximizes likelihood of data given linearized model `lr`.

# Arguments
- `sr::SteadyResults`: Output of [`call_prepare_linearization()`](@ref)
- `lr::LinearResults`: Output of [`linearize_full_model()`](@ref)
- `m_par::ModelParameters`
- `e_set::EstimationSettings`

# Returns
- `EstimResults`, containing all returns of [`mode_finding()`](@ref)
- `posterior_mode::Array{Float64,1}`: Mode of the posterior
- `smoother_output::Array{Float64,2}`: Smoother output
- `sr::SteadyResults`: Updated `SteadyResults`
- `lr::LinearResults`: Updated `LinearResults`
- `m_par::ModelParameters`: Updated `ModelParameters`
"""
function find_mode(
    sr::SteadyResults,
    lr::LinearResults,
    m_par::ModelParameters,
    e_set::EstimationSettings,
)
    @printf "\n"
    @printf "Started mode finding. This might take a while...\n"

    priors = collect(metaflatten(m_par, prior)) # prior distributions of model parameters that are estimated
    if e_set.mode_start_file == ""
        if e_set.me_treatment != :fixed
            append!(priors, e_set.meas_error_distr) # add the meas. error priors
        end
        par_start = mode.(priors)
    else
        # Load the dictionary
        par_final_dict_string = open(e_set.mode_start_file, "r") do file
            JSON.parse(read(file, String)) # Read file and parse JSON
        end
        par_final_dict = Dict(Symbol(k) => v for (k, v) in par_final_dict_string) # Convert string keys back to Symbols

        # Construct the vector par_final dependant on which parameters are estimated
        par_final = Vector{Float64}()
        counter = 1
        for field in fieldnameflatten(m_par)
            value = get(par_final_dict, field, mode(priors[counter])) # Use the mode of m_par prior if intial guess is missing
            par_final = push!(par_final, value)
            counter += 1
        end

        counter = 1
        for field in e_set.meas_error_input
            value = get(par_final_dict, field, mode(e_set.meas_error_distr[counter]))
            par_final = push!(par_final, value)
            counter += 1
        end

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

    @printf "Started mode finding. This might take a while... Done.\n"

    return er, posterior_mode, smoother_output, sr, lr, m_par
end

@doc raw"""
    mcmc_estimation(sr, lr, er, m_par, e_set)

Sample posterior of parameter vector with [`rwmh()`](@ref), take sample mean as parameter
estimate, and save all results in `file`.

# Arguments
- `sr::SteadyResults`: Output of [`call_prepare_linearization()`](@ref)
- `lr::LinearResults`: Output of [`linearize_full_model()`](@ref)
- `er::EstimResults`: Output of [`find_mode()`](@ref)
- `m_par::ModelParameters`
- `e_set::EstimationSettings`

# Returns
- `sr::SteadyResults`: Updated `SteadyResults`
- `lr::LinearResults`: Updated `LinearResults`
- `er::EstimResults`: Updated `EstimResults`
- `m_par::ModelParameters`: Updated `ModelParameters`
- `draws_raw::Array{Float64,2}`: Raw draws from the posterior
- `posterior::Array{Float64,1}`: Posterior
- `accept_rate::Float64`: Acceptance rate
- `par_final::Array{Float64,1}`: Mode of the posterior
- `hessian_sym::Symmetric{Float64,Array{Float64,2}}`: Hessian of the posterior
- `smoother_output::Array{Float64,2}`: Smoother output
"""
function sample_posterior(
    sr::SteadyResults,
    lr::LinearResults,
    er::EstimResults,
    m_par::ModelParameters,
    e_set::EstimationSettings,
)
    @printf "Started MCMC. This might take a while...\n"

    hessian_sym = Symmetric(nearest_spd(inv(er.hessian_final)))
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
            draws_raw[(e_set.burnin + 1):end, :],
            (size(draws_raw[(e_set.burnin + 1):end, :])..., 1),
        ),
        [string(parnames_ascii[i]) for i in eachindex(parnames_ascii)],
    )
    chn_summary = summarize(chn)
    par_final = chn_summary[:, :mean]

    ##
    if e_set.me_treatment != :fixed
        m_par = Flatten.reconstruct(
            m_par,
            par_final[1:(length(par_final) - length(er.meas_error))],
        )
    else
        m_par = Flatten.reconstruct(m_par, par_final)
    end

    lr = update_model(sr, lr, m_par)

    smoother_output = likeli(par_final, sr, lr, er, m_par, e_set; smoother = true)

    @printf "Started MCMC. This might take a while... Done.\n"

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

end
