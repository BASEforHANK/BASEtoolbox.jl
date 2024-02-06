@doc raw"""
    plot_irfs(IRFs, SHOCKs, select_variables, nice_var_names, nice_s_names, horizon, model_names, n_plotcol; savepdf = false, disp_switch = true, suffix = "")

This function is designed to plot impulse response functions (IRFs). 
It takes as input the impulse responses, shocks, selected variables, nice variable names, nice shock names, horizon, model names, and number of plot columns, and generates a plot.

# Arguments
- `IRFs`: The impulse responses to be plotted. This should be a collection of impulse responses (an Array).
- `SHOCKs`: The symbols of the shocks applied to the models.
- `select_variables`: The variables to be displayed in the plot. This should be a collection of variable names, a collection of strings.
- `nice_var_names`: The readable names of the variables to be displayed in the plot. This should be a collection of strings.
- `nice_s_names`: The readable names of the shocks to be displayed in the plot. This should be a collection of strings.
- `horizon`: The horizon over which to plot the impulse response functions.
- `model_names`: The names of the models. This should be an array of strings.
- `n_plotcol`: The number of plot columns.

# Optional Arguments
- `savepdf`: A boolean indicating whether to save the plot as a PDF. Default is false.
- `disp_switch`: A boolean indicating whether to display the switch. Default is true.
- `suffix`: A string to be appended to the end of the file name when saving the plot as a PDF.

# Returns
- This function returns a vector of plots.

# Examples
```julia
plot_irfs(IRFs, SHOCKs, select_variables, nice_var_names, nice_s_names, horizon, model_names, n_plotcol)
```
"""
function plot_irfs(
    IRFs,
    SHOCKs,
    select_variables,
    nice_var_names,
    nice_s_names,
    horizon,
    model_names::Array{String},
    n_plotcol;
    savepdf = false,
    disp_switch = true,
    suffix = "",
)
    styles = [:solid :dash :dashdot :dashdotdot :dot :dash :solid]
    colorlist = [:black, :blue, :firebrick1, :green, :orange, :purple, :yellow]
    pvec = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(SHOCKs))
    counts = 0
    for s in SHOCKs
        counts += 1
        countv = 0
        countm = 0
        p = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(select_variables) + 1)
        for v = 1:length(select_variables)+1
            p[v] = plot()
        end
        lim_l = zeros(length(select_variables) + 1)
        lim_u = zeros(length(select_variables) + 1) .+ 0.001
        for irf in IRFs
            countm += 1
            countv = 1
            p[countv] = plot!(
                p[countv],
                irf[countv, 1:horizon, counts],
                linewidth = 2,
                linestyle = styles[countm],
                palette = colorlist,
                legend = false,
                title = nice_s_names[counts],
                titlefontsize = 10,
                fontfamily = "Computer Modern",
            )
            for v in select_variables
                countv += 1
                lim_l[countv] = min(
                    1.1 * minimum(irf[countv, 1:horizon, counts]) - 0.001,
                    lim_l[countv],
                )
                lim_u[countv] = max(
                    1.1 * maximum(irf[countv, 1:horizon, counts]) + 0.001,
                    lim_u[countv],
                )

                p[countv] = plot!(
                    p[countv],
                    irf[countv, 1:horizon, counts],
                    linewidth = 2,
                    ylims = (lim_l[countv], lim_u[countv]),
                    linestyle = styles[countm],
                    palette = colorlist,
                    legend = false,
                    title = nice_var_names[countv-1],
                    titlefontsize = 10,
                    fontfamily = "Computer Modern",
                )
            end
        end
        # Add a plot that only contains the legend
        p[1] = plot!(
            p[1],
            label = model_names,
            legend = :topright,
            linestyle = styles,
            palette = colorlist,
            legendfontsize = 8,
            foreground_color_legend = nothing,
            background_color_legend = nothing,
            fontfamily = "Computer Modern",
        )
        # Combine in a plot with sublots
        rows = ceil(Int, length(p) / n_plotcol)
        pvec[counts] = plot(
            p...,
            layout = (rows, n_plotcol),
            size = rows .* (400, 250) .* 0.7,
            linewidth = 2,
            thickness_scaling = 1.1,
            fontfamily = "Computer Modern",
            label = model_names,
            tickfontvalign = :top,
        )
    end
    if savepdf
        for j in eachindex(SHOCKs)
            savefig(
                pvec[j],
                string("Output/Figures/IRFs/IRFs_to_", SHOCKs[j], suffix, ".pdf"),
            )
        end
    end
    if disp_switch
        display.(pvec)
    end
    return pvec
end
