"""
    plot_vardecomp(
        vars_to_plot,
        VDs_to_plot,
        VDs_order,
        ids;
        shock_categories = Dict(),
        VD_horizons = [4, 16, 100],
        colorlist = [...],
        show_fig = true,
        save_fig = false,
        path = "",
        suffix = "",
        factor = 100.0,
    )

Plots variance decompositions (VDs) for selected variables across multiple models and
forecast horizons, given a set of VDs as produced by `compute_vardecomp`, potentially
organized by shock categories.

# Arguments

  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, where each tuple
    consists of a variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `VDs_to_plot::Vector{Tuple{Array{Float64,3},String}}`: A vector of tuples, where each
    tuple consists of a 3D array of VDs (`Array{Float64,3}`) and a string representing the
    specification name.
  - `VDs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    VD arrays.
  - `ids`: A structure mapping variable symbols to their corresponding indices in the VD
    arrays, must be identical for all VDs versions.

# Keyword Arguments

  - `shock_categories::Dict{Tuple{String,String},Vector{Symbol}}`: A dictionary where each
    key represents a category of shocks (with a label and a string for saving), and each
    value is a vector of symbols representing the shocks in that category.
  - `VD_horizons::Vector{Int64}`: Forecast horizons to consider. Default is `[4, 16, 100]`.
  - `colorlist::Vector`: A list of colors used to distinguish different shocks in the
    stacked bar plots. Default includes predefined hex colors and symbols.
  - `show_fig::Bool`: If `true`, displays the plot. Default is `true`.
  - `save_fig::Bool`: If `true`, saves the plot as a PDF. Default is `false`.
  - `path::String`: The directory path where the generated plots should be saved. Default is
    an empty string (no saving).
  - `suffix::String`: A suffix to append to the saved plot filenames. Default is an empty
    string.
  - `factor::Float64 = 100.0`: Scaling factor for variance decomposition values.
"""
function plot_vardecomp(
    vars_to_plot::Vector{Tuple{Symbol,String}},
    VDs_to_plot::Vector{Tuple{Array{Float64,3},String}},
    VDs_order::Vector{Symbol},
    ids;
    shock_categories::Dict{Tuple{String,String},Vector{Symbol}} = Dict{
        Tuple{String,String},
        Vector{Symbol},
    }(),
    VD_horizons::Vector{Int64} = [4, 16, 100],
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
    factor::Float64 = 100.0,
)
    # Ensure consistency in data dimensions
    @assert length(colorlist) >= length(VDs_order)

    # Unpack variables and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]
    var_idx = [getfield(ids, vars[i]) for i in eachindex(vars)]

    # Unpack model names
    model_names = [VDs_to_plot[i][2] for i in eachindex(VDs_to_plot)]

    # Check if shock categories are provided, if not, create one category for each shock
    if isempty(shock_categories)
        shock_categories = Dict((string(i), string(i)) => [i] for i in VDs_order)
    end

    # Map shocks to categories
    shock_to_category = Dict(
        shock => category[1] for (category, shocks) in shock_categories for shock in shocks
    )
    ShockCategory = [get(shock_to_category, shock, "Other") for shock in VDs_order]

    # Convert shock and variable names to categorical arrays for plotting
    col_models = CategoricalArray(string.(model_names))
    levels!(col_models, col_models)
    col_shocks = CategoricalArray(string.(VDs_order))
    levels!(col_shocks, col_shocks)
    col_vars = CategoricalArray(string.(vars))
    levels!(col_vars, col_vars)
    col_labs = CategoricalArray(labs)
    levels!(col_labs, col_labs)

    # Get dimensions
    M, H, N, S =
        length(model_names), length(VD_horizons), length(col_vars), length(col_shocks)

    # Convert VDs_to_plot into a DataFrame, select only relevant variables
    VDs_to_plot_DF = vcat(
        [
            DataFrame(;
                Model = model_names[m],
                Horizon = VD_horizons[h],
                Variable = col_vars[n],
                Label = col_labs[n],
                VarianceDecomposition = VDs_to_plot[m][1][var_idx, :, :][
                    n,
                    VD_horizons[h],
                    s,
                ] * factor,
                Shock = ShockCategory[s],
            ) for m = 1:M, h = 1:H, n = 1:N, s = 1:S
        ]...,
    )

    # Aggregate contributions by shock category
    VDs_to_plot_DF = combine(
        groupby(VDs_to_plot_DF, [:Model, :Horizon, :Variable, :Label, :Shock]),
        :VarianceDecomposition => sum,
    )

    # Plot each variable separately
    for (i, i_var) in enumerate(vars)

        # Select data for current variable
        df = VDs_to_plot_DF[VDs_to_plot_DF.Variable .== string(i_var), :]

        # Create stacked bar plot
        p =
            data(df) *
            mapping(
                :Variable,
                :VarianceDecomposition_sum;
                col = :Horizon => nonnumeric,
                row = :Model,
                stack = :Shock,
                color = :Shock,
            ) *
            visual(BarPlot; direction = :x)

        # Draw figure with custom styling
        pp = draw(
            p,
            scales(; Color = (; palette = colorlist));
            axis = (;
                xgridvisible = false,
                ygridvisible = false,
                yticksvisible = false,
                yticklabelsvisible = false,
                limits = (0, 100, nothing, nothing),
                xticks = (0:20:100),
                xlabel = "Shock Contribution in %",
                ylabel = "",
            ),
            figure = (; size = (1600, 400),),
            legend = (; framevisible = false),
        )

        # Add title to figure
        pp.figure[0, 1:end] = Label(
            pp.figure,
            string("Variance Decomposition for ", labs[i], ", Forecast Horizon:");
            fontsize = 20,
        )

        # Save plot
        if save_fig
            AlgebraOfGraphics.save(path * "/CVD_" * string(i_var) * suffix * ".pdf", pp)
        end

        # Show plot
        if show_fig
            display(pp)
        end
    end
end

"""
    plot_vardecomp_bcfreq(
        vars_to_plot,
        VDbcs_to_plot,
        VDbcs_order,
        ids;
        shock_categories = Dict(),
        colorlist = [...],
        show_fig = true,
        save_fig = false,
        path = "",
        suffix = "",
        factor = 100.0,
    )

Plots busincess cycle frequency variance decompositions (VDbcs) for selected variables
across multiple models and forecast horizons, given a set of VDbcs as produced by
`compute_vardecomp_bcfreq`, potentially organized by shock categories.

# Arguments

  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, where each tuple
    consists of a variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `VDbcs_to_plot::Vector{Tuple{Matrix{Float64},String}}`: A vector of tuples, where each
    tuple consists of a matrix of VDbcs (`Matrix{Float64}`) and a string representing the
    specification name.
  - `VDbcs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    VDbc arrays.
  - `ids`: A structure mapping variable symbols to their corresponding indices in the VD
    arrays, must be identical for all VDbcs versions.

# Keyword Arguments

  - `shock_categories::Dict{Tuple{String,String},Vector{Symbol}}`: A dictionary where each
    key represents a category of shocks (with a label and a string for saving), and each
    value is a vector of symbols representing the shocks in that category.
  - `colorlist::Vector`: A list of colors used to distinguish different shocks in the
    stacked bar plots. Default includes predefined hex colors and symbols.
  - `show_fig::Bool`: If `true`, displays the plot. Default is `true`.
  - `save_fig::Bool`: If `true`, saves the plot as a PDF. Default is `false`.
  - `path::String`: The directory path where the generated plots should be saved. Default is
    an empty string (no saving).
  - `suffix::String`: A suffix to append to the saved plot filenames. Default is an empty
    string.
  - `factor::Float64 = 100.0`: Scaling factor for variance decomposition values.
"""
function plot_vardecomp_bcfreq(
    vars_to_plot::Vector{Tuple{Symbol,String}},
    VDbcs_to_plot::Vector{Tuple{Matrix{Float64},String}},
    VDbcs_order::Vector{Symbol},
    ids;
    shock_categories::Dict{Tuple{String,String},Vector{Symbol}} = Dict{
        Tuple{String,String},
        Vector{Symbol},
    }(),
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
    factor::Float64 = 100.0,
)
    # Ensure consistency in data dimensions
    @assert length(colorlist) >= length(VDbcs_order)

    # Unpack variables and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]
    var_idx = [getfield(ids, vars[i]) for i in eachindex(vars)]

    # Unpack model names
    model_names = [VDbcs_to_plot[i][2] for i in eachindex(VDbcs_to_plot)]

    # Check if shock categories are provided, if not, create one category for each shock
    if isempty(shock_categories)
        shock_categories = Dict((string(i), string(i)) => [i] for i in VDbcs_order)
    end

    # Map shocks to categories
    shock_to_category = Dict(
        shock => category[1] for (category, shocks) in shock_categories for shock in shocks
    )
    ShockCategory = [get(shock_to_category, shock, "Other") for shock in VDbcs_order]

    # Convert shock and variable names to categorical arrays for plotting
    col_models = CategoricalArray(string.(model_names))
    levels!(col_models, col_models)
    col_shocks = CategoricalArray(string.(VDbcs_order))
    levels!(col_shocks, col_shocks)
    col_vars = CategoricalArray(string.(vars))
    levels!(col_vars, col_vars)
    col_labs = CategoricalArray(labs)
    levels!(col_labs, col_labs)

    # Get dimensions
    M, N, S = length(model_names), length(col_vars), length(col_shocks)

    # Convert VDbcs_to_plot into a DataFrame, select only relevant variables
    VDbcs_to_plot_DF = vcat(
        [
            DataFrame(;
                Model = model_names[m],
                Variable = col_vars[n],
                Label = col_labs[n],
                VarianceDecomposition = VDbcs_to_plot[m][1][var_idx, :][n, s] * factor,
                Shock = ShockCategory[s],
            ) for m = 1:M, n = 1:N, s = 1:S
        ]...,
    )

    # Aggregate contributions by shock category
    VDbcs_to_plot_DF = combine(
        groupby(VDbcs_to_plot_DF, [:Model, :Variable, :Label, :Shock]),
        :VarianceDecomposition => sum,
    )

    # Plot each variable separately
    for (i, i_var) in enumerate(vars)

        # Select data for current variable
        df = VDbcs_to_plot_DF[VDbcs_to_plot_DF.Variable .== string(i_var), :]

        # Create stacked bar plot
        p =
            data(df) *
            mapping(
                :Variable,
                :VarianceDecomposition_sum;
                row = :Model,
                stack = :Shock,
                color = :Shock,
            ) *
            visual(BarPlot; direction = :x)

        # Draw figure with custom styling
        pp = draw(
            p,
            scales(; Color = (; palette = colorlist));
            axis = (;
                xgridvisible = false,
                ygridvisible = false,
                yticksvisible = false,
                yticklabelsvisible = false,
                limits = (0, 100, nothing, nothing),
                xticks = (0:20:100),
                xlabel = "Shock Contribution in %",
                ylabel = "",
            ),
            figure = (; size = (1600, 400),),
            legend = (; framevisible = false),
        )

        # Add title to figure
        pp.figure[0, 1:end] = Label(
            pp.figure,
            string(
                "Variance Decomposition for ",
                labs[i],
                ", at business cycle frequency:",
            );
            fontsize = 20,
        )

        # Save plot
        if save_fig
            AlgebraOfGraphics.save(path * "/CVDbc_" * string(i_var) * suffix * ".pdf", pp)
        end

        # Show plot
        if show_fig
            display(pp)
        end
    end
end
