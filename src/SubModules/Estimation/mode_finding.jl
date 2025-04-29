@doc raw"""
    mode_finding(XSSaggr, A, B, indexes, indexes_aggr, distrSS, compressionIndexes, m_par, n_par, e_set)

Given definition of observed variables and their transformation (level or growth rate) from `e_set`,
load the data, construct the observation equation, and maximize [`likeli()`](@ref) (the log-likelihood)
using the package `Optim`.

Save estimation results to `e_set.save_mode_file`.

# Returns
- `par_final`: parameter vector that maximizes the likelihood
- `hessian_final`: Hessian of the log-likelihood at `par_final`
- `posterior_mode`: log-likelihood at `par_final`
- `meas_error`,`meas_error_std`: returns from [`measurement_error()`](@ref)
- `parnames`: names of estimated parameters (including measurement error variances)
- `Data`,`Data_missing`: data from `e_set.data_file`; marker for missing data
- `H_sel`: selector matrix for states/controls that are observed
- `priors`: priors of parameters (including measurement error variances)
- `smoother_output`: output from the Kalman smoother
"""
function mode_finding(sr, lr, m_par, e_set, par_start)

    # Load data
    Data_temp = DataFrame(CSV.File(e_set.data_file; missingstring = "NaN"))
    data_names_temp = propertynames(Data_temp)

    # Rename observables that do not have matching model names
    for i in data_names_temp
        name_temp = get(e_set.data_rename, i, :none)
        if name_temp != :none
            rename!(Data_temp, Dict(i => name_temp))
        end
    end

    # Identify missing observations
    observed_vars = e_set.observed_vars_input
    Data = Matrix(Data_temp[:, observed_vars])
    Data_missing = ismissing.(Data)

    # Built selection matrix
    H_sel = zeros(e_set.nobservables, sr.n_par.nstates_r + sr.n_par.ncontrols_r)
    for i in eachindex(observed_vars)
        H_sel[i, getfield(sr.indexes_r, (observed_vars[i]))] = 1.0
    end

    # get names of estimated parameters and add measurement error params
    parnames = collect(fieldnameflatten(m_par))
    if e_set.me_treatment != :fixed
        for i in eachindex(e_set.meas_error_input)
            push!(parnames, Symbol(:σ_me_, e_set.meas_error_input[i]))
        end
    end

    # Set up measurement error
    meas_error, meas_error_prior, meas_error_std =
        measurement_error(Data, observed_vars, e_set)

    # initialize parameters at starting values
    par = copy(par_start)

    if e_set.me_treatment != :fixed
        m_par = Flatten.reconstruct(m_par, par[1:(length(par) - length(meas_error))])
    else
        m_par = Flatten.reconstruct(m_par, par)
    end

    # Prior specification
    priors = collect(metaflatten(m_par, prior)) # model parameters
    if e_set.me_treatment != :fixed
        append!(priors, meas_error_prior)          # add the meas. error priors
    end

    # Optimization
    # Define objective function
    Laux(pp) =
        -likeli(
            pp,
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
        )[3]

    # Code variant with box-constrained optimization, used for updating compression
    OptOpt = Optim.Options(;
        show_trace = true,
        show_every = 20,
        store_trace = true,
        x_abstol = e_set.x_tol,
        f_reltol = e_set.f_tol,
        iterations = div(e_set.max_iter_mode, 3),
    )
    opti = optimize(Laux, par, e_set.optimizer, OptOpt)

    par_final = Optim.minimizer(opti)
    # Update estimated model parameters and resolve model
    if e_set.me_treatment != :fixed
        m_par = Flatten.reconstruct(
            m_par,
            par_final[1:(length(par_final) - length(meas_error))],
        )
    else
        m_par = Flatten.reconstruct(m_par, par_final)
    end
    if sr.n_par.verbose
        @printf "Updating model reduction after initial optimization.\n"
    end
    @set! sr.n_par.further_compress = false
    sr_aux = model_reduction(sr, lr, m_par) # revert to full model
    lr_aux = update_model(sr_aux, lr, m_par) # solve full model
    @set! sr_aux.n_par.further_compress = true
    sr = model_reduction(sr_aux, lr_aux, m_par) # update model reduction
    lr = update_model(sr, lr_aux, m_par)   # solve new reduced model

    # Built selection matrix
    H_sel = zeros(e_set.nobservables, sr.n_par.nstates_r + sr.n_par.ncontrols_r)
    for i in eachindex(observed_vars)
        H_sel[i, getfield(sr.indexes_r, (observed_vars[i]))] = 1.0
    end

    # Redefine objective function
    LL(pp) =
        -likeli(
            pp,
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
        )[3]

    OptOpt = Optim.Options(;
        show_trace = true,
        show_every = 20,
        store_trace = true,
        x_abstol = e_set.x_tol,
        f_reltol = e_set.f_tol,
        iterations = div(e_set.max_iter_mode, 3) * 2,
    )
    opti = optimize(LL, Optim.minimizer(opti), e_set.optimizer, OptOpt)
    par_final = Optim.minimizer(opti)

    # Update estimated model parameters and resolve model
    if e_set.me_treatment != :fixed
        m_par = Flatten.reconstruct(
            m_par,
            par_final[1:(length(par_final) - length(meas_error))],
        )
    else
        m_par = Flatten.reconstruct(m_par, par_final)
    end
    ll_old = -Optim.minimum(opti)

    if sr.n_par.verbose
        @printf "Updating model reduction after mode finding finished.\n"
    end
    @set! sr.n_par.further_compress = false
    sr_aux = model_reduction(sr, lr, m_par) # revert to full model
    lr_aux = update_model(sr_aux, lr, m_par) # solve full model
    if sr.n_par.verbose
        @printf "New reduction.\n"
    end
    @set! sr_aux.n_par.further_compress = true
    sr = model_reduction(sr_aux, lr_aux, m_par) # update model reduction
    lr = update_model(sr, lr_aux, m_par)   # solve new reduced model
    # Built selection matrix
    H_sel = zeros(e_set.nobservables, sr.n_par.nstates_r + sr.n_par.ncontrols_r)
    for i in eachindex(observed_vars)
        H_sel[i, getfield(sr.indexes_r, (observed_vars[i]))] = 1.0
    end

    LL_final(pp) =
        -likeli(
            pp,
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
        )[3]

    posterior_mode = -LL_final(par_final)
    if sr.n_par.verbose
        @printf "Likelihood at mode under reduction: old: %f new: %f\n" ll_old posterior_mode
    end

    # Run Kalman smoother
    smoother_output = likeli(
        par_final,
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
        smoother = true,
    )

    # Compute Hessian at posterior mode
    if e_set.compute_hessian == true
        if sr.n_par.verbose
            @printf "Computing Hessian. This might take a while...\n"
        end
        hessian_final = FiniteDiff.finite_difference_hessian(LL, par_final; relstep = 0.001)
    else
        if sr.n_par.verbose
            @printf "Assuming Hessian is I...\n"
        end
        hessian_final = Matrix{Float64}(I, length(par_final), length(par_final))
    end

    return par_final,
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
    lr
end
