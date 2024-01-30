@doc raw"""
    compute_irfs_vardecomp(models, select_variables)

This function calculates a variance decomposition of the given models for selected variables. 
It returns four outputs: the impulse responses, the fixed horizon decompositions, 
the list of shocks that drive the model, and finally the variance decomposition based on the 
method proposed by Uhlig (2001).

# Arguments
- `models`: The models for which the variance decomposition will be calculated. This should be a collection of models.
- `select_variables`: The variables selected for the variance decomposition. This should be a collection of variable names.

# Returns
- `IRFs`: The impulse responses of the models.
- `VDs`: The fixed horizon decompositions of the models.
- `SHOCKs`: The shocks applied to the models.
- `VD_bc_s`: The variance decomposition based on the method proposed by Harald Uhlig (2001).

# Examples
```julia
IRFs, VDs, SHOCKs, VD_bc_s = compute_irfs_vardecomp(models, select_variables)
```
"""
function compute_irfs_vardecomp(models, select_variables)

    max_horizon = 1000
    n_models = length(models)
    IRFs = Array{Array{Float64}}(undef, n_models)
    VDs = Array{Array{Float64}}(undef, n_models)
    VD_bc_s = Array{Array{Float64}}(undef, n_models)
    SHOCKs = Array{Symbol}(undef, 0)
    for j = 1:n_models
        e_set = models[j][3]
        union!(SHOCKs, e_set.shock_names)
    end
    n_shocks = length(SHOCKs)
    for j = 1:n_models
        sr = models[j][1]
        lr = models[j][2]
        e_set = models[j][3]
        m_par = models[j][4]
        selector = []
        isstate = zeros(Bool, length(select_variables) + 1)
        iter = 2
        for i in select_variables
            if i in Symbol.(sr.state_names)
                isstate[iter] = true
            end
            iter += 1
            try
                append!(selector, getfield(sr.indexes_r, i))
            catch
                append!(selector, sr.n_par.ntotal_r + 1)
                println(selector)
            end
        end

        IRFs_aux = zeros(length(selector) + 1, max_horizon, n_shocks)
        shock_number = 0
        for i in SHOCKs
            x = zeros(size(lr.LOMstate, 1))
            shock_number += 1
            selectorplus = copy(selector)
            try
                shock = getfield(sr.indexes_r, i)
                x[shock] = getfield(m_par, Symbol("σ_", i))
                selectorplus = [shock; selector]
            catch
                println("model ", j, " has no shock ", string(i))
            end
            MX = [I; lr.State2Control; zeros(1, sr.n_par.nstates_r)]
            for t = 1:max_horizon
                IRFs_aux[:, t, shock_number] = MX[selectorplus, :] * x
                x[:] = lr.LOMstate * x
            end
        end
        #IRFs_aux[isstate, 1:end-1, :] .= IRFs_aux[isstate, 2:end, :] # IRFs for state variables represent end-of-period values
        # Dimensions of IRF: variable x time x shock
        # VARdecomp = zeros(size(IRFs_aux))

        VARdecomp =
            cumsum(IRFs_aux[2:end, :, :] .^ 2.0, dims = 2) ./
            (sum(cumsum(IRFs_aux[2:end, :, :] .^ 2.0, dims = 2), dims = 3) .+ 10.0 * eps())

        VDs[j] = 100.0 .* VARdecomp
        IRFs[j] = 100.0 .* IRFs_aux
        VD_bc_aux, _ =
            compute_bcfreq_vardecomp(sr, lr, e_set, m_par; passband = (6, 32), ngrid = 512)
        VD_bc_aux = [VD_bc_aux; zeros(1, size(VD_bc_aux, 2))]
        VD_bc_s[j] = VD_bc_aux[selector, :] .* 100.0
    end
    return IRFs, VDs, SHOCKs, VD_bc_s

end
