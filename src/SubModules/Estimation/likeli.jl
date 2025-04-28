@doc raw"""
    likeli(par, Data, Data_missing, H_sel, priors, meas_error, meas_error_std, sr, lr, m_par, e_set; smoother=false)
Compute the likelihood of `Data`, given model-parameters `par` and prior `priors` (maximize to find MLE of `par`).

Solve model with [`LinearSolution_reduced_system()`](@ref), compute likelihood with [`kalman_filter()`](@ref) or with [`kalman_filter_smoother()`](@ref) (if `smoother==True`).

# Returns
*if `smoother==False`:*
- `log_like`,`prior_like`,`post_like`,`alarm`: log-likelihoods (`post` is the sum of `prior` and computed likelihood); `alarm` indicates error when solving model with [`LinearSolution_reduced_system`](@ref), sets log-likelihood to `-9.e15`
*if `smoother==True`:*
- `smoother_output`: returns from [`kalman_filter_smoother()`](@ref)
"""
function likeli(
    par,
    Data,
    Data_missing,
    H_sel,
    priors,
    meas_error,
    meas_error_std,
    sr,
    lr,
    m_par,
    e_set;
    smoother = false,
)
    return likeli_backend(
        par,
        Data,
        Data_missing,
        H_sel,
        priors,
        meas_error,
        meas_error_std,
        sr,
        lr,
        m_par,
        e_set,
        smoother,
    )
end

@doc raw"""
    likeli(par, sr, lr, er, m_par, e_set; smoother=false)
Compute the likelihood of `er.Data`, given model-parameters `par` and prior `er.priors` (maximize to find MLE of `par`).

Solve model with [`LinearSolution_reduced_system()`](@ref), compute likelihood with [`kalman_filter()`](@ref) or with [`kalman_filter_smoother()`](@ref) (if `smoother==True`).

# Returns
*if `smoother==False`:*
- `log_like`,`prior_like`,`post_like`,`alarm`: log-likelihoods (`post` is the sum of `prior` and computed likelihood); `alarm` indicates error when solving model with [`LinearSolution_reduced_system`](@ref), sets log-likelihood to `-9.e15`
*if `smoother==True`:*
- `smoother_output`: returns from [`kalman_filter_smoother()`](@ref)
- `State2Control`,`LOM`: state-to-control and state transition matrizzes
"""
function likeli(par, sr, lr, er, m_par, e_set; smoother = false)
    return likeli_backend(
        par,
        er.Data,
        er.Data_missing,
        er.H_sel,
        er.priors,
        er.meas_error,
        er.meas_error_std,
        sr,
        lr,
        m_par,
        e_set,
        smoother,
    )
end

function likeli_backend(
    par,
    Data,
    Data_missing,
    H_sel,
    priors,
    meas_error,
    meas_error_std,
    sr,
    lr,
    m_par,
    e_set,
    smoother,
)

    # check priors, abort if they are violated
    prior_like::eltype(par), alarm_prior::Bool = prioreval(Tuple(par), Tuple(priors))
    alarm = false
    if alarm_prior
        log_like = -9.e15
        alarm = true
        State2Control = zeros(sr.n_par.ncontrols_r, sr.n_par.nstates_r)
        if e_set.debug_print
            @printf "Parameter try violates PRIOR.\n"
        end
    else
        if e_set.me_treatment != :fixed
            m_start = length(par) - length(meas_error) # find out where in par structural pars end
        else
            m_start = length(par)
        end

        # replace estimated values in m_par by last candidate
        m_par = Flatten.reconstruct(m_par, par[1:m_start])

        # covariance of structural shocks
        SCov = zeros(eltype(par), sr.n_par.nstates_r, sr.n_par.nstates_r)
        for i in e_set.shock_names
            SCov[getfield(sr.indexes_r, i), getfield(sr.indexes_r, i)] =
                (getfield(m_par, Symbol("σ_", i))) .^ 2
        end

        # covariance of measurement errors, assumption: ME ordered after everything else
        m = size(H_sel)[1]
        MCov = diagm(zeros(eltype(par), m)) # no correlated ME allowed for now
        if !isempty(meas_error)
            m_iter = 1
            if e_set.me_treatment != :fixed
                for (k, v) in meas_error # read out position of measurement errors
                    MCov[v, v] = par[m_start + m_iter] .^ 2
                    m_iter += 1
                end
            else
                for (k, v) in meas_error # read out position of measurement errors
                    MCov[v, v] = meas_error_std[m_iter] .^ 2
                    m_iter += 1
                end
            end
        end

        # solve model using candidate parameters
        # BLAS.set_num_threads(1)
        State2Control::Array{eltype(par),2},
        LOM::Array{eltype(par),2},
        alarm_LinearSolution::Bool =
            LinearSolution_reduced_system(sr, m_par, lr.A, lr.B; allow_approx_sol = false)

        # BLAS.set_num_threads(Threads.nthreads())
        if alarm_LinearSolution # abort if model doesn't solve
            log_like = -9.e15
            alarm = true
            if e_set.debug_print
                @printf "Parameter try leads to inexistent or unstable equilibrium.\n"
            end
        else
            MX = [I; State2Control]
            H = H_sel * MX
            if smoother == false
                log_like = kalman_filter(H, LOM, Data, Data_missing, SCov, MCov, e_set)
                # log_like = kalman_filter_herbst(Data, LOM, SCov, H, MCov, 0, e_set)
            else
                smoother_output =
                    kalman_filter_smoother(H, LOM, Data, .!Data_missing, SCov, MCov, e_set)
                log_like = smoother_output[1]
            end
        end
    end
    post_like = log_like + prior_like

    if smoother == false
        return log_like, prior_like, post_like, alarm, State2Control
    else
        return smoother_output
    end
end

@doc raw"""
    prioreval(par,priors)

Evaluate prior PDF at the parameters given in `par`.

# Arguments
- `par`: vector of parameters [npar*1]
- `priors`: vector of prior distributions [npar*1]

# Returns
- `log_priorval`: log prior density [scalar]
- `alarm`: indicator that is 1 if there is a violation of the prior bounds [scalar]
"""
function prioreval(par, priors)
    if all(insupport.(priors, par))
        alarm = false
        log_priorval = sum(logpdf.(priors, par))
    else
        alarm = true
        log_priorval = -9.e15
    end

    return log_priorval, alarm
end
