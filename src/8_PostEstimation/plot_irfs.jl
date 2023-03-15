###############################################################################################
# Compute IRFs and variance decomposition for set of models and variables passed to function
###############################################################################################




###############################################################################################
# Plot IRFs
###############################################################################################
function plot_irfs(IRFs, SHOCKs, select_variables, nice_var_names, nice_s_names, horizon, model_names::Array{String}, n_plotcol; savepdf = false, suffix="")
    styles = [:solid :dash :dashdot :dashdotdot :dot :dash :solid]
    colorlist = [ :black, :blue, :firebrick1, :green, :orange, :purple, :yellow]
    pvec = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(SHOCKs))
    counts = 0
    for s in SHOCKs
        counts += 1
        countv = 0
        countm = 0
        p = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(select_variables) + 1)
        for v=1:length(select_variables)+1
            p[v] = plot()
        end
        lim_l=zeros(length(select_variables) + 1)
        lim_u=zeros(length(select_variables) + 1).+.001
        for irf in IRFs
            countm += 1
            countv = 1
            p[countv] = plot!(p[countv], 
                            irf[countv, 1:horizon, counts], 
                            linewidth = 2,
                            linestyle = styles[countm], 
                            palette = colorlist, 
                            legend = false, title = nice_s_names[counts],
                            titlefontsize=10,  
                            fontfamily = "Computer Modern"
                            )
            for v in select_variables
                countv += 1       
                lim_l[countv] = min(1.1*minimum(irf[countv, 1:horizon, counts]) -.001, lim_l[countv])
                lim_u[countv] = max(1.1*maximum(irf[countv, 1:horizon, counts]) +.001, lim_u[countv])
            
                p[countv] = plot!(p[countv], 
                                    irf[countv, 1:horizon, counts], 
                                    linewidth = 2,
                                    ylims =(lim_l[countv],lim_u[countv]),
                                    linestyle = styles[countm], 
                                    palette = colorlist, 
                                    legend = false, title = nice_var_names[countv-1],
                                    titlefontsize=10, 
                                    fontfamily = "Computer Modern"
                                    )
            end
        end
        # Add a plot that only contains the legend
        p[1] = plot!(p[1], 
                        label = model_names, 
                        legend = :topright,
                        linestyle = styles, 
                        palette = colorlist, 
                        legendfontsize = 8, 
                        foreground_color_legend = nothing,
                        background_color_legend = nothing, 
                        fontfamily = "Computer Modern")
        # Combine in a plot with sublots
        rows =ceil(Int, length(p) / n_plotcol)
        pvec[counts] = plot(p..., 
                        layout = (rows, n_plotcol), 
                        size = rows.*(400, 250).*0.7,
                        linewidth = 2, 
                        thickness_scaling = 1.1, 
                        fontfamily = "Computer Modern", 
                        label = model_names,
                        tickfontvalign=:top,
                        )
    end
    if savepdf
        for j = eachindex(SHOCKs)
            savefig(pvec[j],string("8_PostEstimation/Figures/IRFs/IRFs_to_", SHOCKs[j],suffix, ".pdf"))
        end
    end
    display.(pvec)
    return pvec
end
