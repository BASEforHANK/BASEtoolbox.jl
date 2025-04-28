"""
    plot_irfs_cat(
        shock_categories,
        vars_to_plot,
        IRFs_to_plot,
        IRFs_order,
        ids;
        horizon = 40,
        factor = 100,
        show_fig = true,
        save_fig = false,
        save_fig_indiv = false,
        path = "",
        suffix = "",
        yscale = "standard",
        style_options = (lw = 2, color = :auto, linestyle = :solid)
    )

Plots impulse response functions (IRFs) for specified shocks and variables, given IRFs as
produced by `compute_irfs`, organized by shock categories.

# Arguments

  - `shock_categories::Dict{Tuple{String,String},Vector{Symbol}}`: A dictionary where each
    key represents a category of shocks (with a label and a string for saving), and each
    value is a vector of symbols representing the shocks in that category.
  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, each containing a
    variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `IRFs_to_plot::Array{Float64,3}`: A 3D array of IRFs.
  - `IRFs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    IRF arrays.
  - `ids`: A structure mapping variable symbols to their corresponding indices in the IRF
    arrays, which must be identical for all IRFs versions.

# Keyword Arguments

  - `horizon::Int64`: The time horizon (number of periods) over which IRFs are plotted.
    Default is `40`.
  - `factor::Int64`: Scaling factor for the IRFs (default: `100`).
  - `show_fig::Bool`: If `true`, displays the plot. Default is `true`.
  - `save_fig::Bool`: If `true`, saves the combined plot as a PDF. Default is `false`.
  - `save_fig_indiv::Bool`: If `true`, saves individual plots for each variable or shock
    panel directly inside a category-specific folder. Default is `false`.
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
function plot_irfs_cat(
    shock_categories::Dict{Tuple{String,String},Vector{Symbol}},
    vars_to_plot::Vector{Tuple{Symbol,String}},
    IRFs_to_plot::Array{Float64,3},
    IRFs_order::Vector{Symbol},
    ids;
    horizon::Int64 = 40,
    factor::Int64 = 100,
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
        lw = 2,
        dpi = 300,
        size = (1600, 1000),
        foreground_color_legend = nothing,
        background_color_legend = nothing,
    )

    # Unpack variables and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]

    # Define the base directory for the IRFs folder
    irfs_path = joinpath(path)
    mkpath(irfs_path)

    # Loop over categories
    for (category_name, category_shocks) in shock_categories

        # Find position of current shocks (category_shocks) in IRFs array (IRFs_order)
        idx = [findfirst(x -> x == i_shock, IRFs_order) for i_shock in category_shocks]

        # If one of the shocks is not found, skip the category and print a warning
        if any(isnothing, idx)
            @warn "One or more shocks in category $category_name not found in IRFs_order, skipped category."
            continue
        end

        # Extract IRFs for these shocks, round to 10 digits, and multiply by 100
        i_IRFs = mapround(IRFs_to_plot[:, :, idx]; digits = 10) .* factor
        n_IRFs = size(i_IRFs, 3)

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

        # Determine y-axis limits for "common" yscale, tuple, or dictionary
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
        for (i, (var, lab)) in enumerate(zip(vars, labs))
            if hasfield(typeof(ids), var)
                var_idx = getfield(ids, var)
                p = plot(
                    i_IRFs[var_idx, 1:horizon, 1];
                    title = lab,
                    label = string(category_shocks[1]),
                    lw = style_options.lw,
                    color = effective_color[1],
                    linestyle = effective_linestyle[1],
                    tickfont = font(10, "Computer Modern"),
                    titlefont = font(12, "Computer Modern"),
                )
                for j = 2:n_IRFs
                    plot!(
                        p,
                        i_IRFs[var_idx, 1:horizon, j];
                        label = string(category_shocks[j]),
                        lw = style_options.lw,
                        color = effective_color[j],
                        linestyle = effective_linestyle[j],
                    )
                end
                # Apply y-axis limits based on dictionary values
                if yscale isa Dict && haskey(yscale, var)
                    ylims!(p, ylimits_per_variable[var]...)
                elseif yscale == "common" || yscale isa Tuple{Number,Number}
                    ylims!(p, ymin, ymax)
                end
                p = plot!(p; legend = false)
            else
                @printf "Variable %s not found in ids\n" var
                p = plot()
            end
            push!(pp, p)
        end

        # Add panel for the shocks
        p = plot(;
            title = category_name[1],
            lw = style_options.lw,
            linestyle = style_options.linestyle,
            tickfont = font(10, "Computer Modern"),
            titlefont = font(12, "Computer Modern"),
        )
        for (i, var) in enumerate(category_shocks)
            if hasfield(typeof(ids), var)
                var_idx = getfield(ids, var)
                plot!(
                    p,
                    i_IRFs[var_idx, 1:horizon, i];
                    label = string(var),
                    lw = style_options.lw,
                    color = effective_color[i],
                    linestyle = effective_linestyle[i],
                )
            else
                @printf "Variable %s not found in ids\n" var
            end
        end
        pp = [p; pp...]

        # Combine all plots in a single figure, with a layout of nrow x ncol
        n = length(pp)
        ncol = ceil(Int, sqrt(n))
        nrow = ceil(Int, n / ncol)
        fig = plot(pp...; layout = (nrow, ncol), pp_layout...)

        # Save combined plot directly in the IRFs folder
        if save_fig
            savefig(fig, joinpath(irfs_path, "IRF_" * category_name[2] * suffix * ".pdf"))
        end

        # Save individual plots directly inside each category folder
        if save_fig_indiv
            category_path = joinpath(irfs_path, category_name[2])
            mkpath(category_path)
            for (i, p) in enumerate(pp)
                p = plot!(p; legend = true, pp_layout...)
                var_name = i > length(vars) ? "ShockPanel" : string(vars[i])  # Handle extra subplot
                savefig(
                    p,
                    joinpath(
                        category_path,
                        "IRF_" * category_name[2] * "_" * var_name * suffix * ".pdf",
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
