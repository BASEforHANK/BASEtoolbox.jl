"""
    plot_hist_decomp(
        vars_to_plot,
        HDs_to_plot,
        HDs_order,
        ids;
        shock_categories = Dict(),
        timeline = collect(1:size(HDs_to_plot, 2)),
        colorlist = [...],
        show_fig = true,
        save_fig = false,
        path = "",
        suffix = ""
    )

Plots historical decompositions (HDs) for specified variables, showing the contribution of
different shocks over time, given a set of HDs as produced by `compute_hist_decomp`,
potentially organized by shock categories.

# Arguments

  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, where each tuple
    consists of a variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `HDs_to_plot::Array{Float64,3}`: A 3D array containing historical decompositions, where
    dimensions correspond to variables, time periods, and shocks.
  - `HDs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    HDs array.
  - `ids`: A structure mapping variable symbols to their corresponding indices in the HDs
    array, ensuring consistency across all plotted variables.

# Keyword Arguments

  - `shock_categories::Dict{Tuple{String,String},Vector{Symbol}}`: A dictionary where each
    key represents a category of shocks (with a label and a string for saving), and each
    value is a vector of symbols representing the shocks in that category.
  - `timeline::Vector`: A vector specifying the time axis for the plots. Default is a
    sequence from `1` to the number of time periods in `HDs_to_plot`.
  - `colorlist::Vector`: A list of colors used to distinguish different shocks in the
    stacked bar plots. Default includes predefined hex colors and symbols.
  - `show_fig::Bool`: If `true`, displays the plot. Default is `true`.
  - `save_fig::Bool`: If `true`, saves the plot as a PDF. Default is `false`.
  - `path::String`: The directory path where the generated plots should be saved. Default is
    an empty string (no saving).
  - `suffix::String`: A suffix to append to the saved plot filenames. Default is an empty
    string.
"""
function plot_hist_decomp(
    vars_to_plot::Vector{Tuple{Symbol,String}},
    HDs_to_plot::Array{Float64,3},
    HDs_order::Vector{Symbol},
    ids;
    shock_categories::Dict{Tuple{String,String},Vector{Symbol}} = Dict{
        Tuple{String,String},
        Vector{Symbol},
    }(),
    timeline::Vector = collect(1:size(HDs_to_plot, 2)),
    colorlist::Vector = [
        "#ff0000",
        "#ffcccc",
        "#ff9900",
        "#ffebcc",
        "#009900",
        "#ccffcc",
        :indigo,
        "#0000ff",
        "#ccccff",
        :grey,
    ],
    show_fig::Bool = true,
    save_fig::Bool = false,
    path::String = "",
    suffix::String = "",
)

    # Ensure consistency in data dimensions
    @assert size(HDs_to_plot, 2) == length(timeline)
    @assert length(colorlist) >= length(HDs_order)

    # Unpack variables and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]
    var_idx = [getfield(ids, vars[i]) for i in eachindex(vars)]

    # Check if shock categories are provided, if not, create one category for each shock
    if isempty(shock_categories)
        shock_categories = Dict((string(i), string(i)) => [i] for i in HDs_order)
    end

    # Map shocks to categories
    shock_to_category = Dict(
        shock => category[1] for (category, shocks) in shock_categories for shock in shocks
    )
    shock_to_category[:initial] = "Initial"
    ShockCategory = [get(shock_to_category, shock, "Other") for shock in HDs_order]

    # Convert shock and variable names to categorical arrays for plotting
    col_shocks = CategoricalArray(string.(HDs_order))
    levels!(col_shocks, col_shocks)
    col_vars = CategoricalArray(string.(vars))
    levels!(col_vars, col_vars)
    col_labs = CategoricalArray(labs)
    levels!(col_labs, col_labs)

    # Get dimensions
    T, N, S = length(timeline), length(col_vars), length(col_shocks)

    # Convert HDs_to_plot into a DataFrame, select only relevant variables
    HDs_to_plot_DF = vcat(
        [
            DataFrame(;
                Time = timeline[t],
                Variable = col_vars[n],
                Label = col_labs[n],
                Contribution = HDs_to_plot[var_idx, :, :][n, t, s],
                Shock = ShockCategory[s],
            ) for t = 1:T, n = 1:N, s = 1:S
        ]...,
    )

    # Aggregate contributions by shock category
    HDs_to_plot_DF = combine(
        groupby(HDs_to_plot_DF, [:Time, :Variable, :Label, :Shock]),
        :Contribution => sum,
    )

    # Plot each variable separately
    for (i, i_var) in enumerate(vars)

        # Select data for current variable
        df = HDs_to_plot_DF[HDs_to_plot_DF.Variable .== string(i_var), :]
        df_tot = combine(groupby(df, :Time), :Contribution_sum => sum)

        # Stacked bar plot of shock contributions over time
        pp = @df df groupedbar(
            :Time,
            :Contribution_sum,
            group = :Shock,
            bar_position = :stack,
            linecolor = false,
            palette = colorlist,
            size = (400, 300),
        )

        # Overlay total contribution as a black line
        pp = @df df_tot plot!(
            :Time,
            :Contribution_sum_sum,
            linewidth = 1,
            linecolor = :black,
            label = "total",
        )

        # Adjust x-axis limits
        xlims!((minimum(timeline) - 2, maximum(timeline) + 2))

        # Add title to plot
        plot!(pp; title = string(col_labs[i]), fontfamily = "Computer Modern")

        # Add legend to plot
        plot!(
            pp;
            legend = :best,
            legend_columns = 3,
            legendfontsize = 8,
            foreground_color_legend = nothing,
            background_color_legend = nothing,
        )

        # Save plot
        if save_fig
            savefig(pp, path * "/HistD_" * string(i_var) * suffix * ".pdf")
        end

        # Show plot
        if show_fig
            display(pp)
        end
    end
end
