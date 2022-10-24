###############################################################################################
# Compute IRFs and variance decomposition for set of models and variables passed to function
###############################################################################################


function compute_vardecomp_bounds(models, select_variables, model_names; n_replic = 1000, percentile_bounds = (0.05, 0.95))

    _,_, SHOCKs, VDorig = compute_irfs_vardecomp(models, select_variables)

    n_models = length(models)
    VD_lower = Array{Array{Float64}}(undef, n_models)
    VD_upper = Array{Array{Float64}}(undef, n_models)

    n_shocks = length(SHOCKs)
    for j = 1:n_models
        println(j)
        sr = models[j][1]
        lr = models[j][2]
        e_set = models[j][3]
        m_par = models[j][4]
        draws = models[j][5]
        selector = []
        isstate = zeros(Bool, length(select_variables))
        iter = 1
        for i in select_variables
            if i in Symbol.(sr.state_names)
                isstate[iter] = true
            end
            iter += 1
            try
                append!(selector, getfield(sr.indexes_r, i))
            catch
                append!(selector, sr.ntotal_r + 1)
            end
        end
        VDaux = hcat([zeros(n_replic) for j = 1:length(select_variables), k = 1:n_shocks])
        draw_ind = rand(1:size(draws, 1), n_replic)
        Threads.@threads for s = eachindex(draw_ind)
            A=lr.A; B=lr.B; State2Control=lr.State2Control; LOMstate=lr.LOMstate; SolutionError=lr.SolutionError; nk=lr.nk
            lr_local = LinearResults(copy(State2Control), copy(LOMstate), copy(A), copy(B), copy(SolutionError), copy(nk))
            dd = draws[draw_ind[s], :]
            if e_set.me_treatment != :fixed
                m_par_local = Flatten.reconstruct(m_par, dd[1:size(draws, 2)-length(e_set.meas_error_input)])
            else
                m_par_local = Flatten.reconstruct(m_par, dd)
            end
            lr_local = update_model(sr, lr_local, m_par_local)
            VD_bc_s, _ = compute_bcfreq_vardecomp(sr, lr_local, e_set, m_par_local; passband=(6,32), ngrid=512)
            # Dimensions of IRF: variable x time x shock
    
            for d1 = 1:length(select_variables)
                for d2 = 1:n_shocks
                    VDaux[d1, d2][s] = VD_bc_s[selector[d1], d2]*100
                end
            end
        end
    
        VD_lower[j] = quantile.(VDaux, percentile_bounds[1])
        VD_upper[j] = quantile.(VDaux, percentile_bounds[2])
    end
    n_total_entries = n_models * length(select_variables) * n_shocks
    modelname_vec = Vector{String}(undef, n_total_entries)
    variablename_vec = Vector{Symbol}(undef, n_total_entries)
    shockname_vec = Vector{Symbol}(undef, n_total_entries)
    VD_vec = Vector{Float64}(undef, n_total_entries)
    VDlow_vec = Vector{Float64}(undef, n_total_entries)
    VDup_vec = Vector{Float64}(undef, n_total_entries)
    count = 1
    for d0 = 1:n_models
        for d1 = 1:length(select_variables)
            for d2 = 1:n_shocks
                modelname_vec[count] = model_names[d0]
                variablename_vec[count] = select_variables[d1]
                shockname_vec[count] = SHOCKs[d2]
                VD_vec[count] = copy(VDorig[d0][d1, d2])
                VDlow_vec[count] = copy(VD_lower[d0][d1, d2])
                VDup_vec[count] = copy(VD_upper[d0][d1, d2])
                count += 1
            end
        end
    end
    VDdf = DataFrame(model = modelname_vec, variable = variablename_vec, shock = shockname_vec, 
                     lower_bound = VDlow_vec, point_estimate = VD_vec, upper_bound = VDup_vec)

    return VDdf

end
