@doc raw"""
    plot_vardecomp(VDs, VD_bc_s, select_vd_horizons, model_names, SHOCKs, select_variables; savepdf = false, suffix = "", legend_switch = true, disp_switch = true)

This function is designed to plot variance decompositions. It takes as input the variance decompositions. 
These are the variance decomposition based on the method proposed by Harald Uhlig (2001) as well as variance 
decomposition at fixed horizons. It further expects as inputs model names, shocks as symbols, and the variables
selected to be displayed. It generates a plot.

# Arguments
- `VDs`: The variance decompositions to be plotted. This should be a collection of variance decompositions at fixed horizons.
- `VD_bc_s`: The variance decomposition based on the method proposed by Harald Uhlig (2001).
- `select_vd_horizons`: The horizons of the variance decompositions to be displayed in the plot. This should be a collection of integers.
- `model_names`: The names of the models for which the varaiance decompositions are provided. This should be an array of strings.
- `SHOCKs`: The symbols of the shocks applied to the models.
- `select_variables`: The variables to be displayed in the plot. This should be a collection of variable names.

# Optional Arguments
- `savepdf`: A boolean indicating whether to save the plot as a PDF. Default is false.
- `suffix`: A string to be appended to the end of the file name when saving the plot as a PDF.
- `legend_switch`: A boolean indicating whether to display the legend. Default is true.
- `disp_switch`: A boolean indicating whether to display the switch. Default is true.

# Returns
- This function returns the plotted variance decompositions as data frames.

# Examples
```julia
plot_vardecomp(VDs, VD_bc_s, select_vd_horizons, model_names, SHOCKs, select_variables)
```
"""
function plot_vardecomp(
    VDs,
    VD_bc_s,
    select_vd_horizons,
    model_names,
    SHOCKs,
    select_variables;
    savepdf = false,
    suffix = "",
    legend_switch = true,
    disp_switch = true
)
    shock_color = [
        "#ff0000",
        "#ffcccc",
        "#ff9900",
        "#ffebcc",
        "#009900",
        "#ccffcc",
        :indigo,
        "#0000ff",
        "#ccccff",
    ]
    SHOCKs = CategoricalArray(string.(SHOCKs))
    levels!(SHOCKs, SHOCKs)
    select_variables = CategoricalArray(string.(select_variables))
    levels!(select_variables, select_variables)
    model_names = CategoricalArray(model_names[:])
    levels!(model_names, model_names)
    # construct DataFrame that stacks variance decomposition
    df = vcat(
        [
            DataFrame(
                Model = model_names[i],
                Horizon = h,
                VarianceDecomposition = VDs[i][j, h, k],
                Shock = SHOCKs[k],
                Variable = select_variables[j],
            ) for i in eachindex(model_names), h in select_vd_horizons,
            k in eachindex(SHOCKs), j in eachindex(select_variables)
        ]...,
    )

    df_bc = vcat(
        [
            DataFrame(
                Model = model_names[i],
                VarianceDecomposition = VD_bc_s[i][j, k],
                Shock = SHOCKs[k],
                Variable = select_variables[j],
            ) for i in eachindex(model_names), k in eachindex(SHOCKs),
            j in eachindex(select_variables)
        ]...,
    )

    # plot variance decomposition
    for j in select_variables
        p = ((
            df[df.Variable.==j, :] |> @vlplot(
                :bar,
                x = {
                    :VarianceDecomposition,
                    title = nothing,
                    scale = {domain = [0, 100], nice = false},
                },
                color = {
                    :Shock,
                    scale = {range = shock_color},
                    sort = SHOCKs,
                    legend = legend_switch,
                },
                row = {"Model:n", title = "", sort = model_names},
                column = {
                    "Horizon:n",
                    title = string("Variance Decomposition for ", j, ", Forecast Horizon:"),
                    sort = select_vd_horizons,
                },
                background = :white,
                order = "siteOrder:o",
                width = 220,
                height = 50
            )
        ))
        if savepdf
            save(string("Output/Figures/CVD/CVD_of_", j, suffix, ".pdf"), p)
        end
        if disp_switch 
            display(p)
        end
    end
    for j in select_variables
        p = ((
            df_bc[df_bc.Variable.==j, :] |> @vlplot(
                :bar,
                x = {
                    :VarianceDecomposition,
                    title = nothing,
                    scale = {domain = [0, 100], nice = false},
                },
                color = {
                    :Shock,
                    scale = {range = shock_color},
                    sort = SHOCKs,
                    legend = legend_switch,
                },
                row = {"Model:n", title = "", sort = model_names},
                background = :white,
                order = "siteOrder:o",
                width = 220,
                height = 55
            )
        ))
        if savepdf
            save(
                string("Output/Figures/CVD/CVD_at_bcf_of_", j, suffix, ".pdf"),
                p,
            )
        end
        if disp_switch
            display(p)
        end
    end

    return df, df_bc
end
