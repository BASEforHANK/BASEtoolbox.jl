"""
    mode_finding(sr, lr, m_par, e_set, par_start)

Find the posterior mode of the parameter vector via numerical optimization.

Two estimation modes are supported:

  - Likelihood estimation (default, when `e_set.irf_matching == false`): loads time-series
    data, builds the observation matrix, and maximizes the log posterior via `likeli` using
    `Optim`.
  - IRF matching (when `e_set.irf_matching == true`): loads IRF targets from CSV according
    to `e_set.irf_matching_dict`, builds weights, and maximizes the IRF-matching objective
    via `irfmatch`.

In both modes, the routine updates the reduced model using `model_reduction`/`update_model`
between optimization rounds to improve accuracy, and it optionally computes a Hessian at the
mode.

# Arguments

  - `sr`: structural/reduced-system container (sizes, indexes, reduction flags)
  - `lr`: linearization container (A, B matrices and caches)
  - `m_par`: model parameter struct (Flatten-compatible)
  - `e_set::EstimationSettings`: controls data paths, observables, priors, optimizer and
    tolerances, measurement-error treatment, IRF-matching settings, and flags like
    `compute_hessian`
  - `par_start::AbstractVector`: initial parameter vector (includes measurement-error stds
    if estimated)

# Returns

Returns a 16-tuple with elements depending on the mode:

 1. `par_final::Vector{Float64}`: parameter vector at the mode
 2. `hessian_final::Matrix{Float64}`: Hessian at the mode (finite-diff or identity if
    disabled)
 3. `posterior_mode::Float64`: value of the objective at the mode (log posterior)
 4. `meas_error`: mapping from observable names to column indexes (likelihood mode); `[]`
    for IRF matching
 5. `meas_error_std::Vector{Float64}`: std caps or fixed stds (likelihood mode); `[]` for
    IRF matching
 6. `parnames::Vector{Symbol}`: names of estimated parameters (incl. ME params if not fixed)
 7. `Data::Matrix{Float64}`: data matrix (likelihood mode); `[]` for IRF matching
 8. `Data_missing::BitMatrix`: missingness mask (likelihood mode); `[]` for IRF matching
 9. `IRFtargets::Array{Float64,3}`: target IRFs (IRF matching); empty otherwise
10. `IRFserrors::Array{Float64,3}`: IRF target standard errors (IRF matching); empty
    otherwise
11. `H_sel::Matrix{Float64}`: selection matrix mapping states/controls to observables
    (likelihood mode); empty for IRF matching
12. `priors::Vector{<:Distribution}`: priors for structural (and ME) parameters
13. `smoother_output`: Kalman smoother output from `likeli(...; smoother=true)` (likelihood
    mode); empty for IRF matching
14. `m_par`: parameter struct updated at `par_final`
15. `sr`: structural container after final reduction update
16. `lr`: linearization container after final update

Notes:

  - Measurement errors are constructed by `measurement_error` and included in `parnames` and
    `priors` unless `e_set.me_treatment == :fixed`. # No irf matching? => likelihood
    estimation
  - The optimizer and tolerances are taken from `e_set` (e.g., `e_set.optimizer`,
    `e_set.x_tol`).        # Load data
"""
function mode_finding(sr, lr, m_par, e_set, par_start; skip_optimization::Bool = false)
    if !e_set.irf_matching # No irf matching? => likelihood estimation
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

        if skip_optimization
            @printf "skip_optimization=true: skipping numerical optimization, returning prior modes.\n"
            par_final = mode.(priors)
            hessian_final = Matrix{Float64}(I, length(par_final), length(par_final))

            return par_final,
            hessian_final,
            NaN,
            meas_error,
            meas_error_std,
            parnames,
            Data,
            Data_missing,
            [],
            [],
            H_sel,
            priors,
            [],
            m_par,
            sr,
            lr
        end

        # Optimization
        @printf "Running numerical optimization...\n"
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
            hessian_final =
                FiniteDiff.finite_difference_hessian(LL_final, par_final; relstep = 0.001)
        else
            if sr.n_par.verbose
                @printf "Assuming Hessian is I...\n"
            end
            hessian_final = Matrix{Float64}(I, length(par_final), length(par_final))
        end

        # For IRF matching, but we do not use it here
        IRFtargets = []
        IRFserrors = []

        return par_final,
        hessian_final,
        posterior_mode,
        meas_error,
        meas_error_std,
        parnames,
        Data,
        Data_missing,
        IRFtargets,
        IRFserrors,
        H_sel,
        priors,
        smoother_output,
        m_par,
        sr,
        lr

    elseif e_set.irf_matching # irf matching? => irf matching estimation
        if sr.n_par.verbose
            @printf "IRF matching...\n"
        end

        # initialize parameters at starting values
        par = copy(par_start) #[1:(end - length(e_set.meas_error_distr))]

        parnames = collect(fieldnameflatten(m_par))

        # Prior specification
        priors = collect(metaflatten(m_par, prior)) # model parameters

        irf_horizon = e_set.irf_matching_dict["irf_horizon"]

        # Load irf targets
        key_to_use =
            haskey(e_set.irf_matching_dict, "irfs_to_target") ? "irfs_to_target" :
            "irfs_to_plot"

        Data_temp =
            DataFrame(CSV.File(e_set.irf_matching_dict[key_to_use]; missingstring = "NaN"))
        shock_names = Symbol.(unique(Data_temp[:, :shock]))

        data_names_temp = Symbol.(names(Data_temp))

        # Rename observables that do not have matching model names
        for i in data_names_temp
            name_temp = get(e_set.data_rename, i, :none)
            if name_temp != :none
                rename!(Data_temp, Dict(i => name_temp))
            end
        end

        shocks_selected = intersect(shock_names, e_set.shock_names)
        select_variables =
            intersect(Symbol.(propertynames(Data_temp)), e_set.observed_vars_input)

        IRFtargets = Array{Float64}(
            undef,
            irf_horizon,
            length(select_variables),
            length(shocks_selected),
        )
        IRFserrors = Array{Float64}(
            undef,
            irf_horizon,
            length(select_variables),
            length(shocks_selected),
        )

        counti = 0
        for i in shocks_selected
            counti += 1
            countj = 0
            for j in select_variables
                countj += 1
                IRFtargets[:, countj, counti] = Data_temp[
                    (Data_temp[
                        :,
                        :pointdum,
                    ] .== 1) .& (Symbol.(Data_temp[:, :shock]) .== i),
                    j,
                ]
                IRFserrors[:, countj, counti] = Data_temp[
                    (Data_temp[
                        :,
                        :pointdum,
                    ] .== 0) .& (Symbol.(Data_temp[:, :shock]) .== i),
                    j,
                ]
            end
        end
        var_to_scale_by = e_set.irf_matching_dict["scale_responses_by"]
        scale_term =
            var_to_scale_by === nothing ? 1 : maximum(Data_temp[:, var_to_scale_by])

        IRFtargets = IRFtargets ./ scale_term
        IRFserrors = IRFserrors ./ scale_term

        weights = 1.0 ./ (IRFserrors .^ 2)
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

        if skip_optimization
            @printf "skip_optimization=true: skipping numerical optimization, returning prior modes.\n"
            par_final = mode.(priors)
            hessian_final = Matrix{Float64}(I, length(par_final), length(par_final))

            return par_final,
            hessian_final,
            NaN,
            [],
            [],
            parnames,
            [],
            [],
            IRFtargets,
            IRFserrors,
            [],
            priors,
            [],
            m_par,
            sr,
            lr
        end

        m_par = Flatten.reconstruct(m_par, par)

        # Optimization
        @printf "Running numerical optimization (IRF matching)...\n"
        # Define objective function
        Lirfaux(pp) =
            -irfmatch(
                pp,
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
            )[3]

        OptOpt = Optim.Options(;
            show_trace = true,
            show_every = 20,
            store_trace = true,
            x_abstol = e_set.x_tol,
            f_reltol = e_set.f_tol,
            iterations = e_set.max_iter_mode,
        )

        opti = optimize(
            Lirfaux,
            par_start, # [1:(end - length(e_set.meas_error_distr))],
            e_set.optimizer,
            OptOpt,
        )
        par_final = Optim.minimizer(opti)

        m_par = Flatten.reconstruct(m_par, par_final)
        # ll_old = -Optim.minimum(opti)

        if sr.n_par.verbose
            @printf "Updating model reduction after initial optimization.\n"
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

        LLirf_final(pp) =
            -irfmatch(
                pp,
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
            )[3]

        posterior_mode = -LLirf_final(par_final)
        # if sr.n_par.verbose
        #     @printf "Likelihood at mode under reduction: old: %f new: %f\n" ll_old posterior_mode
        # end

        # Compute Hessian at posterior mode
        if e_set.compute_hessian == true
            if sr.n_par.verbose
                @printf "Computing Hessian. This might take a while...\n"
            end
            func = TwiceDifferentiable(pp -> LLirf_final(pp), par_final)
            hessian_final = Optim.hessian!(func, par_final)
        else
            if sr.n_par.verbose
                @printf "Assuming Hessian is I...\n"
            end
            hessian_final = Matrix{Float64}(I, length(par_final), length(par_final))
        end

        meas_error = []
        meas_error_std = []
        H_sel = []
        smoother_output = []
        Data = []
        Data_missing = []

        return par_final,
        hessian_final,
        posterior_mode,
        meas_error,
        meas_error_std,
        parnames,
        Data,
        Data_missing,
        IRFtargets,
        IRFserrors,
        H_sel,
        priors,
        smoother_output,
        m_par,
        sr,
        lr

    else
        error("estimation type not defined")
    end
end
