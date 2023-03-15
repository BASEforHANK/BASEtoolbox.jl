###############################################################################################
# Compute IRFs and variance decomposition for set of models and variables passed to function
###############################################################################################
###############################################################################################
# Plot historical decomposition
###############################################################################################
function compute_hist_decomp(sr, lr, e_set, m_par, smoother_output, select_variables,timeline; 
                            savepdf=false, prefix="")
    SHOCKs = e_set.shock_names
    n_shocks = length(SHOCKs)
    n_select_var = length(select_variables)
    T = size(smoother_output[6], 2)
    IRFs = Array{Float64}(undef, n_select_var, T, n_shocks + 1)
    
    # select variables of interest from all variables
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
            append!(selector, sr.n_par.ntotal_r + 1)
        end
    end

    MX = [I; lr.State2Control]

    # effect of initial condition
    IRFs_aux = zeros(length(selector), T)
    x = smoother_output[3][:,1] 
    for t = 1:T
        IRFs_aux[:, t] = MX[selector, :] * x
        x[:] = lr.LOMstate * x
    end
    IRFs[:,:,n_shocks+1] = IRFs_aux

    # loop through shocks, calculate effect of each shock
    for j  = 1:n_shocks
        IRFs_aux = zeros(length(selector), T)
        x = zeros(size(lr.LOMstate, 1))
        i = SHOCKs[j] 
        shock_index = getfield(sr.indexes_r, i)
        for t = 1:T
            IRFs_aux[:, t] = MX[selector, :] * x
            x[:] = lr.LOMstate * x
            x[shock_index] += smoother_output[6][shock_index,t] # shock in "t" moves observables in "t+1" 
        end
        IRFs[:,:,j] = IRFs_aux
    end

    SHOCKsPlus = push!(copy(SHOCKs), :initial)
    SHOCKsPlus = CategoricalArray(string.(SHOCKsPlus))
    levels!(SHOCKsPlus, SHOCKsPlus)
    select_variables = CategoricalArray(string.(select_variables))
    levels!(select_variables ,select_variables )
    HistDecDF  = vcat([DataFrame(Time = timeline[t], 
                                Variable = select_variables[v], 
                                Contribution = IRFs[v,t,s], 
                                Shock = SHOCKsPlus[s]) for t=1:T, v=1:n_select_var, s=1:n_shocks+1]...)
    
    colorlist = ["#ff0000", "#ffcccc", "#ff9900", "#ffebcc", "#009900", "#ccffcc", :indigo, "#0000ff", "#ccccff", :grey]

    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(select_variables) + 1)
    # rec = vspan(recessions_vec, fill=:grey85, linecolor=:grey85, label="")
    countj=0
    for j in select_variables
        countj+=1

        p[countj] = @df HistDecDF[HistDecDF.Variable .== j,:] groupedbar(:Time, :Contribution, group = :Shock, linecolor =false, title = j,bar_position = :stack, legend=false, palette = colorlist, size=(400,300))
        p[countj] = @df combine(groupby(HistDecDF[HistDecDF.Variable .== j,:], :Time), :Contribution => sum) plot!(:Time, :Contribution_sum, label = string(j), linewidth =1, linecolor=:black, legend=false, fontfamily = "Computer Modern")
        Plots.xlims!((minimum(timeline)-2, maximum(timeline)+2))
    end
    p[end] = bar(Matrix{Missing}(undef, 3, n_shocks+1), label = [string(z) for k = 1:1, z in SHOCKsPlus],linecolor =false, legend = :inside,
                    palette = colorlist, framestyle = :none, legendfontsize = 10)
    p[end] = plot!(Matrix{Missing}(undef, 3, 1), label = "Total", legend = :inside,
                    linecolor=:black, framestyle = :none, legendfontsize = 10, fontfamily = "Computer Modern")
    if savepdf
        # plot(p...) |> save(string("8_PostEstimation/Figures/Smooth/HistD.pdf"))
        for j = eachindex(select_variables)
            plot!(p[j], title="", fontfamily = "Computer Modern") 
            savefig(p[j],string("8_PostEstimation/Figures/Smooth/HistD_",prefix, select_variables[j], ".pdf"))
            display(p[j])

        end
    end
    display(plot(p...))
    return IRFs, HistDecDF, p
end

