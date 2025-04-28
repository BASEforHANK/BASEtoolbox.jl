"""
    plot_irfs(
        shocks_to_plot,
        vars_to_plot,
        IRFs_to_plot,
        IRFs_order,
        ids;
        horizon = 40,
        factor = 100,
        legend_on_all = false,
        show_fig = true,
        save_fig = false,
        save_fig_indiv = false,
        path = "",
        suffix = "",
        yscale = "standard",
        style_options = (lw = 2, color = :auto, linestyle = :solid),
    )

Plots impulse response functions (IRFs) for specified shocks and variables, given a set of
IRFs as produced by `compute_irfs`.

# Arguments

  - `shocks_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, where each tuple
    contains a shock variable (as a `Symbol`) and its corresponding label (`String`).
  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, each containing a
    variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `IRFs_to_plot::Vector{Tuple{Array{Float64,3},String}}`: A vector of tuples, where each
    tuple consists of a 3D array of IRFs (`Array{Float64,3}`) and a string representing the
    specification name.
  - `IRFs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    IRF arrays.
  - `ids`: A structure mapping variable symbols to their corresponding indices in the IRF
    arrays, must be identical for all IRFs versions.

# Keyword Arguments

  - `horizon::Int64`: The time horizon (number of periods) over which IRFs are plotted.
    Default is `40`.
  - `factor::Int64`: Scaling factor for the IRFs (default: `100`).
  - `legend_on_all::Bool`: If `true`, includes the legend on all subplots; otherwise, only
    the first subplot has a legend. Default is `false`.
  - `show_fig::Bool`: If `true`, displays the plot. Default is `true`.
  - `save_fig::Bool`: If `true`, saves the combined plot as a PDF. Default is `false`.
  - `save_fig_indiv::Bool`: If `true`, saves individual plots for each variable directly
    inside a shock-specific folder. Default is `false`.
  - `path::String`: The directory path where the generated plots should be saved. Default is
    an empty string (no saving).
  - `suffix::String`: A suffix to append to the saved plot filenames. Default is an empty
    string.
  - `yscale::Union{String, Tuple{Number,Number}, Dict{Symbol,Tuple{Number,Number}}}`: Y-axis
    scaling specification. When set to `"common"`, computes a common y-axis limit across all
    subplots from the data; when provided as a tuple, uses it as `(ymin, ymax)`; when
    provided as a dictionary, applies specified y-axis limits for each variable. Default is
    `"standard"`, which applies default scaling.
  - `style_options::NamedTuple`: A named tuple specifying stylistic options for the plots,
    including line width (`lw`), color (default: `:auto`), and linestyle (default:
    `:solid`). Default is `(lw = 2, color = :auto, linestyle = :solid)`.
"""
function plot_irfs(
    shocks_to_plot::Vector{Tuple{Symbol,String}},
    vars_to_plot::Vector{Tuple{Symbol,String}},
    IRFs_to_plot::Vector{Tuple{Array{Float64,3},String}},
    IRFs_order::Vector{Symbol},
    ids;
    horizon::Int64 = 40,
    factor::Int64 = 100,
    legend_on_all::Bool = false,
    show_fig::Bool = true,
    save_fig::Bool = false,
    save_fig_indiv::Bool = false,
    path::String = "",
    suffix::String = "",
    yscale::Union{String,Tuple{Number,Number},Dict{Symbol,Tuple{Number,Number}}} = "standard",
    style_options::NamedTuple = (lw = 2, color = :auto, linestyle = :solid),
)

    # General stylistic choices for the plots
    pp_layout = (
        dpi = 300,
        size = (1600, 1000),
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        tickfont = font(10, "Computer Modern"),
        titlefont = font(12, "Computer Modern"),
        lw = style_options.lw,
    )

    # Get number of different IRF specifications (1 for each parameter vector)
    n_IRFs = length(IRFs_to_plot)

    effective_color = if style_options.color == :auto
        pcol = palette(:auto)
        length(pcol) < n_IRFs ? vcat(pcol, fill(pcol[1], n_IRFs - length(pcol))) :
        pcol[1:n_IRFs]
    elseif isa(style_options.color, AbstractVector)
        length(style_options.color) < n_IRFs ?
        vcat(
            style_options.color,
            fill(style_options.color[1], n_IRFs - length(style_options.color)),
        ) : style_options.color[1:n_IRFs]
    else
        fill(style_options.color, n_IRFs)
    end

    effective_linestyle = if isa(style_options.linestyle, AbstractVector)
        length(style_options.linestyle) < n_IRFs ?
        vcat(
            style_options.linestyle,
            fill(style_options.linestyle[1], n_IRFs - length(style_options.linestyle)),
        ) : style_options.linestyle[1:n_IRFs]
    else
        fill(style_options.linestyle, n_IRFs)
    end

    # Get size of IRFs array
    ntotal, T, _ = size(IRFs_to_plot[1][1])

    # Unpack variables (fields) and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]

    # Define the base directory for the IRFs folder
    irfs_path = joinpath(path)
    mkpath(irfs_path)

    for (i_shock, i_shock_lab) in shocks_to_plot

        # Find position of current shock (i_shock) in IRFs array (IRFs_order)
        idx = findfirst(x -> x == i_shock, IRFs_order)

        # If the shock is not found, skip and print a warning
        if isnothing(idx)
            @warn "The shock $i_shock not found in IRFs_order, skipped shock."
            continue
        end

        # Extract IRFs for this shock, round to 10 digits, and multiply by factor
        i_IRFs = fill(NaN64, ntotal, T, n_IRFs)
        for i = 1:n_IRFs
            i_IRFs[:, :, i] = mapround(IRFs_to_plot[i][1][:, :, idx]; digits = 10) .* factor
        end

        # Extract the variables to plot and their labels
        i_vars = [i_shock, vars...]
        i_labs = [i_shock_lab, labs...]

        # Determine y-axis limits for "common" yscale or user-defined tuple
        if yscale == "common"
            ymin, ymax = extrema(
                vcat(
                    [
                        i_IRFs[getfield(ids, var), 1:horizon, :] for
                        var in vars if hasfield(typeof(ids), var)
                    ]...,
                ),
            )
        elseif yscale isa Tuple{Number,Number}
            ymin, ymax = yscale
        elseif yscale isa Dict
            ylimits_per_variable = yscale
        else
            ymin, ymax = nothing, nothing
        end

        # Create a plot for each variable
        pp = []
        for (i, (var, lab)) in enumerate(zip(i_vars, i_labs))
            if hasfield(typeof(ids), var)
                var_idx = getfield(ids, var)
                p = plot(
                    i_IRFs[var_idx, 1:horizon, 1];
                    title = lab,
                    label = IRFs_to_plot[1][2],
                    lw = pp_layout.lw,
                    color = effective_color[1],
                    linestyle = effective_linestyle[1],
                )
                for j = 2:n_IRFs # Code to handle color vector for multiple IRFs
                    plot!(
                        p,
                        i_IRFs[var_idx, 1:horizon, j];
                        label = IRFs_to_plot[j][2],
                        lw = pp_layout.lw,
                        color = effective_color[j],
                        linestyle = effective_linestyle[j],
                    )
                end
                if i > 1 && !legend_on_all
                    p = plot!(p; legend = false)
                end
                if yscale isa Dict && haskey(yscale, var)
                    ylims!(p, ylimits_per_variable[var]...)
                elseif yscale == "common" || yscale isa Tuple{Number,Number}
                    ylims!(p, ymin, ymax)
                end
            else
                @printf "Variable %s not found in ids\n" var
                p = plot()
            end
            push!(pp, p)
        end

        # Combine all plots in a single figure, with a layout of nrow x ncol
        n = length(pp)
        ncol = ceil(Int, sqrt(n))
        nrow = ceil(Int, n / ncol)
        fig = plot(pp...; layout = (nrow, ncol), pp_layout...)

        # Save combined plot directly in the IRFs folder
        if save_fig
            savefig(fig, joinpath(irfs_path, "IRF_" * string(i_shock) * suffix * ".pdf"))
        end

        # Save individual plots directly inside each shock folder
        if save_fig_indiv
            shock_path = joinpath(irfs_path, string(i_shock))
            mkpath(shock_path)
            for (i, p) in enumerate(pp)
                p = plot!(p; legend = true, pp_layout...)
                savefig(
                    p,
                    joinpath(
                        shock_path,
                        "IRF_" *
                        string(i_shock) *
                        "_" *
                        string(i_vars[i]) *
                        suffix *
                        ".pdf",
                    ),
                )
            end
        end

        # Show plot
        if show_fig
            display(fig)
        end
    end
end
