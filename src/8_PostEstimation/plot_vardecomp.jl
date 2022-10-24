###############################################################################################
# Compute IRFs and variance decomposition for set of models and variables passed to function
###############################################################################################

###############################################################################################
# Plot variance decomposition
###############################################################################################
function plot_vardecomp(VDs, VD_bc_s, select_vd_horizons, model_names, SHOCKs, select_variables; savepdf = false, suffix="", legend_switch=true)
    shock_color = ["#ff0000", "#ffcccc", "#ff9900", "#ffebcc", "#009900", "#ccffcc", :indigo, "#0000ff", "#ccccff"]
    SHOCKs = CategoricalArray(string.(SHOCKs))
    levels!(SHOCKs,SHOCKs)
    select_variables = CategoricalArray(string.(select_variables))
    levels!(select_variables,select_variables)
    model_names = CategoricalArray(model_names[:])
    levels!(model_names,model_names)
    # construct DataFrame that stacks variance decomposition
    df = vcat([DataFrame(Model = model_names[i], Horizon = h,
        VarianceDecomposition = VDs[i][j, h, k], Shock = SHOCKs[k], Variable = select_variables[j]) 
        for i = 1:length(model_names), h in select_vd_horizons, k = 1:length(SHOCKs), j = 1:length(select_variables)]...)

    df_bc = vcat([DataFrame(Model = model_names[i],  VarianceDecomposition = VD_bc_s[i][j, k], Shock = SHOCKs[k], Variable = select_variables[j]) 
        for i = 1:length(model_names), k = 1:length(SHOCKs), j = 1:length(select_variables)]...)

    # plot variance decomposition
    for j in select_variables
        p = ((df[df.Variable .== j, :] |> @vlplot(:bar,
            x = {:VarianceDecomposition,
                title = nothing,
                scale = {domain = [0, 100], nice = false}},
            color = {:Shock, scale = {range = shock_color}, sort = SHOCKs, legend = legend_switch},
            row = {"Model:n", title="", sort = model_names},
            column = {"Horizon:n", title = string("Variance Decomposition for ", j, ", Forecast Horizon:"), sort = select_vd_horizons},
            background = :white,
            order = "siteOrder:o",
            width = 220,
            height = 50)))
        if savepdf
            p |> save(string("8_PostEstimation/Figures/CVD/CVD_of_", j, suffix,".pdf"))
        end
        display(p)
    end
    for j in select_variables
        p = ((df_bc[df_bc.Variable .== j, :] |> @vlplot(:bar,
            x = {:VarianceDecomposition,
                title = nothing,
                scale = {domain = [0, 100], nice = false}},
            color = {:Shock, scale = {range = shock_color}, 
                     sort = SHOCKs, legend = legend_switch},
            row = {"Model:n", title = "", sort = model_names},
            background = :white,
            order = "siteOrder:o",
            width = 220,
            height = 35)))
        if savepdf
            p |> save(string("8_PostEstimation/Figures/CVD/CVD_at_bcf_of_", j,suffix, ".pdf"))
        end
        display(p)
    end

    return df
end
