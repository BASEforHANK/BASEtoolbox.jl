"""
    dime_mcmc(xhat, Σ, sr, lr, er, m_par, e_set)

DIME (Differential-Independence Mixture Ensemble) MCMC sampler for posterior inference.

Wraps `DIMESampler.RunDIME` to sample from the posterior of the model parameters, as an
alternative to `rwmh`. DIME runs an ensemble of chains with adaptive proposals and does not
require manual tuning of the proposal covariance.

# Arguments

  - `xhat::Vector{Float64}`: initial parameter vector (posterior mode)
  - `Σ::Symmetric{Float64,Array{Float64,2}}`: covariance estimate (inverse Hessian at mode),
    used to initialize the ensemble
  - `sr, lr, er, m_par, e_set`: model structures and estimation settings

# Returns

  - `draws::Matrix{Float64}`: `(ndraws * nchain) × length(xhat)` matrix of post-burnin draws
  - `posterior::Vector{Float64}`: log posterior for each draw
  - `accept_rate::Float64`: estimated acceptance rate across all chains
"""
function dime_mcmc(
    xhat::Vector{Float64},
    Σ::Symmetric{Float64,Array{Float64,2}},
    sr,
    lr,
    er,
    m_par,
    e_set,
)
    ndim = length(xhat)
    nchain = e_set.dime_nchain_factor * ndim
    niter = e_set.ndraws + e_set.burnin

    # Suppress per-evaluation debug output — too noisy with nchain parallel evaluations
    @set! e_set.debug_print = false

    if !e_set.irf_matching
        vec_logprob = _make_vec_logprob_likeli(sr, lr, er, m_par, e_set)
    else
        IRFtargets, weights, shocks_selected, isstate, indexes_sel_vars, priors =
            _setup_irfmatch_data(sr, m_par, e_set)
        vec_logprob = _make_vec_logprob_irfmatch(
            IRFtargets,
            weights,
            shocks_selected,
            isstate,
            indexes_sel_vars,
            priors,
            sr,
            lr,
            m_par,
            e_set,
        )
    end

    orig_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)

    init = _dime_ensemble_init(xhat, Σ, nchain, vec_logprob)

    @printf "DIME ensemble initialized: %d chains, %d iterations (%d burnin + %d production), %d threads\n" nchain niter e_set.burnin e_set.ndraws Threads.nthreads()
    flush(stdout)
    flush(stderr)

    try
        chains_3d, lprobs_2d, _ = RunDIME(
            vec_logprob,
            init,
            niter;
            sigma = e_set.dime_sigma,
            gamma = e_set.dime_gamma,
            aimh_prob = e_set.dime_aimh_prob,
            df_proposal_dist = e_set.dime_df_proposal_dist,
            rho = e_set.dime_rho,
            progress = true,
        )

        # Discard burnin iterations
        chains_post = chains_3d[(e_set.burnin + 1):end, :, :]   # [ndraws, nchain, ndim]
        lprobs_post = lprobs_2d[(e_set.burnin + 1):end, :]       # [ndraws, nchain]

        npost = size(chains_post, 1)
        draws = reshape(chains_post, npost * nchain, ndim)
        posterior = reshape(lprobs_post, npost * nchain)
        accept_rate = _estimate_accept_rate(chains_3d)

        return draws, posterior, accept_rate
    finally
        LinearAlgebra.BLAS.set_num_threads(orig_blas_threads)
    end
end

function _make_vec_logprob_likeli(sr, lr, er, m_par, e_set)
    function vec_logprob(params_matrix::AbstractMatrix)
        nchain = size(params_matrix, 2)
        lprobs = Vector{Float64}(undef, nchain)
        Threads.@threads for j = 1:nchain
            par_j = @view params_matrix[:, j]
            result = likeli(par_j, sr, lr, er, m_par, e_set)
            post_like, alarm = result[3], result[4]
            lprobs[j] = alarm ? -1e10 : post_like
        end
        return lprobs
    end
    return vec_logprob
end

function _make_vec_logprob_irfmatch(
    IRFtargets,
    weights,
    shocks_selected,
    isstate,
    indexes_sel_vars,
    priors,
    sr,
    lr,
    m_par,
    e_set,
)
    function vec_logprob(params_matrix::AbstractMatrix)
        nchain = size(params_matrix, 2)
        lprobs = Vector{Float64}(undef, nchain)
        Threads.@threads for j = 1:nchain
            par_j = copy(params_matrix[:, j])
            result = irfmatch(
                par_j,
                IRFtargets,
                weights,
                shocks_selected,
                isstate,
                indexes_sel_vars,
                priors,
                sr,
                lr,
                m_par,
                e_set,
            )
            post_like, alarm = result[3], result[4]
            lprobs[j] = alarm ? -1e10 : post_like
        end
        return lprobs
    end
    return vec_logprob
end

function _setup_irfmatch_data(sr, m_par, e_set)
    irf_horizon = e_set.irf_matching_dict["irf_horizon"]

    key_to_use =
        haskey(e_set.irf_matching_dict, "irfs_to_target") ? "irfs_to_target" :
        "irfs_to_plot"

    Data_temp =
        DataFrame(CSV.File(e_set.irf_matching_dict[key_to_use]; missingstring = "NaN"))
    shock_names_local = Symbol.(unique(Data_temp[:, :shock]))

    data_names_temp = Symbol.(names(Data_temp))
    for i in data_names_temp
        name_temp = get(e_set.data_rename, i, :none)
        if name_temp != :none
            rename!(Data_temp, Dict(i => name_temp))
        end
    end

    shocks_selected = intersect(shock_names_local, e_set.shock_names)
    select_variables =
        intersect(Symbol.(propertynames(Data_temp)), e_set.observed_vars_input)

    IRFtargets = Array{Float64}(
        undef,
        irf_horizon,
        length(select_variables),
        length(shocks_selected),
    )
    IRF_sderr = Array{Float64}(
        undef,
        irf_horizon,
        length(select_variables),
        length(shocks_selected),
    )

    count_shock = 0
    for i in shocks_selected
        count_shock += 1
        count_outcm = 0
        for j in select_variables
            count_outcm += 1
            IRFtargets[:, count_outcm, count_shock] = Data_temp[
                (Data_temp[:, :pointdum] .== 1) .& (Symbol.(Data_temp[:, :shock]) .== i),
                j,
            ]
            IRF_sderr[:, count_outcm, count_shock] = Data_temp[
                (Data_temp[:, :pointdum] .== 0) .& (Symbol.(Data_temp[:, :shock]) .== i),
                j,
            ]
        end
    end

    var_to_scale_by = e_set.irf_matching_dict["scale_responses_by"]
    scale_term = var_to_scale_by === nothing ? 1 : maximum(Data_temp[:, var_to_scale_by])

    IRFtargets = IRFtargets ./ scale_term
    IRF_sderr = IRF_sderr ./ scale_term

    weights = 1.0 ./ (IRF_sderr .^ 2)

    iter = 1
    indexes_sel_vars = []
    isstate = zeros(Bool, length(select_variables))
    for i in select_variables
        if i in Symbol.(sr.state_names)
            isstate[iter] = true
        end
        iter += 1
        append!(indexes_sel_vars, getfield(sr.indexes_r, i))
    end

    priors = collect(metaflatten(m_par, prior))

    return IRFtargets, weights, shocks_selected, isstate, indexes_sel_vars, priors
end

function _dime_ensemble_init(
    xhat::Vector{Float64},
    Σ::Symmetric{Float64,Array{Float64,2}},
    nchain::Int,
    vec_logprob::Function;
    max_attempts::Int = 10,
    init_scale::Float64 = 0.5,
)
    ndim = length(xhat)
    NormDist = MvNormal(zeros(ndim), Matrix(Σ))

    init = xhat .+ init_scale .* rand(NormDist, nchain)
    lprobs = vec_logprob(init)

    for attempt = 1:max_attempts
        bad_idx = findall(lprobs .< -1e6)
        if isempty(bad_idx)
            return init
        end
        scale = init_scale / (2^attempt)
        init[:, bad_idx] .= xhat .+ scale .* rand(NormDist, length(bad_idx))
        lprobs[bad_idx] .= vec_logprob(@view init[:, bad_idx])
    end

    # last resort: place remaining bad chains very close to the mode
    bad_idx = findall(lprobs .< -1e6)
    if !isempty(bad_idx)
        for j in bad_idx
            init[:, j] .= xhat .+ 1e-4 .* randn(ndim)
        end
        lprobs[bad_idx] .= vec_logprob(@view init[:, bad_idx])

        still_bad = findall(lprobs .< -1e6)
        if !isempty(still_bad)
            @printf(
                "Warning: %d of %d initial ensemble members have log-prob below -1e6.\n",
                length(still_bad),
                nchain,
            )
        end
    end

    return init
end

function _estimate_accept_rate(chains_3d::Array{Float64,3})
    niter, nchain, _ = size(chains_3d)
    moved = 0
    total = 0
    for j = 1:nchain
        for i = 2:niter
            total += 1
            if @view(chains_3d[i, j, :]) != @view(chains_3d[i - 1, j, :])
                moved += 1
            end
        end
    end
    return moved / total
end
