"""
    plot_distributional_irfs(
        shocks_to_plot,
        vars_to_plot,
        IRFs_to_plot,
        IRFs_order;
        horizon = 40,
        factor = 100,
        legend = false,
        show_fig = true,
        save_fig = false,
        path = "",
        suffix = ""
    )

Plots impulse response functions (IRFs) for specified shocks and distributional variables,
given a set of IRFs as produced by `compute_irfs`.

# Arguments

  - `shocks_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, where each tuple
    contains a shock variable (as a `Symbol`) and its corresponding label (`String`).

  - `vars_to_plot::Vector{Tuple{Symbol,String}}`: A vector of tuples, each containing a
    variable to plot (as a `Symbol`) and its corresponding label (`String`).
  - `IRFs_to_plot::Dict{String,Array{Float64}}`: A dictionary of distributional IRFs, where
    the keys are variable names and the values are arrays of IRFs.
  - `IRFs_order::Vector{Symbol}`: A vector of symbols specifying the order of shocks in the
    IRF arrays.

# Keyword Arguments

  - `horizon::Int64`: The time horizon (number of periods) over which IRFs are plotted.
    Default is `40`.
  - `factor::Int64`: Scaling factor for the IRFs (default: `100`).
  - `legend_on_all::Bool`: If `true`, includes the legend on all subplots. Otherwise, only
  - `show_fig::Bool`: If `true`, displays the plot (default: `true`). the first subplot has
    a legend (default: `false`).
  - `save_fig::Bool`: If `true`, saves the plot as a PDF. Default is `false`.
  - `path::String`: The directory path where the generated plots should be saved. Default is
    an empty string (no saving).
  - `suffix::String`: A suffix to append to the saved plot filenames. Default is an empty
    string.
"""
function plot_distributional_irfs(
    shocks_to_plot::Vector{Tuple{Symbol,String}},
    vars_to_plot::Vector{Tuple{String,String}},
    IRFs_to_plot::Dict{String,Array{Float64}},
    IRFs_order::Vector{Symbol};
    horizon::Int64 = 40,
    legend::Bool = false,
    show_fig::Bool = true,
    save_fig::Bool = false,
    path::String = "",
    suffix::String = "",
)

    # Unpack variables (fields) and labels
    vars = [vars_to_plot[i][1] for i in eachindex(vars_to_plot)]
    labs = [vars_to_plot[i][2] for i in eachindex(vars_to_plot)]

    # Create plots for each shock
    for (i_shock, i_shock_lab) in shocks_to_plot

        # Find position of current shock (i_shock) in IRFs array (IRFs_order)
        idx = findfirst(x -> x == i_shock, IRFs_order)

        # If the shock is not found, skip and print a warning
        if isnothing(idx)
            @warn "The shock $i_shock not found in IRFs_order, skipped shock."
            continue
        end

        # Create a plot for each variable
        for (i, (var, lab)) in enumerate(zip(vars, labs))
            if var in keys(IRFs_to_plot)
                if var != "distr"
                    @printf "Univariate plot for %s\n" var

                    # Select IRFs for the variable, shock, and horizon
                    i_IRFs = IRFs_to_plot[var]
                    i_IRFs = i_IRFs[:, 1:horizon, idx]

                    # Create plot
                    p = plot_smooth_ridge(i_IRFs; legend = legend)
                    plot!(p; title = lab * " for shock to " * i_shock_lab)

                    # Save plot
                    if save_fig
                        savefig(
                            p,
                            path *
                            "/DistIRFs_" *
                            string(i_shock) *
                            "_" *
                            var *
                            suffix *
                            ".pdf",
                        )
                    end

                    # Show plot
                    if show_fig
                        display(p)
                    end

                elseif var == "distr"
                    @printf "Bivariate plot for %s\n" var

                    # Select IRFs for the variable, shock, and horizon
                    i_IRFs = IRFs_to_plot[var]
                    i_IRFs = i_IRFs[:, :, 1:horizon, idx]

                    # Create plot
                    anim = plot_bivariate_animation(
                        IRFs_to_plot[var][:, :, 1:horizon, idx];
                        legend = legend,
                    )

                    # Save plot
                    if save_fig
                        anim = gif(
                            anim,
                            path *
                            "/DistIRFs_" *
                            string(i_shock) *
                            "_" *
                            var *
                            suffix *
                            ".gif";
                            fps = 2,
                            show_msg = false,
                        )
                    else
                        anim = gif(anim; fps = 2, show_msg = false)
                    end

                    # Show plot
                    if show_fig
                        display(anim)
                    end
                else
                    @printf "Variable %s not found in IRFs\n" var
                end
            else
                @printf "Variable %s not found in IRFs\n" var
            end
        end
    end
end

"""
    plot_smooth_ridge(data; n_grid_points = 200, legend = false)

Generates a smooth ridge plot of impulse response function (IRF) distributions over time.

The function first computes kernel density estimates (KDEs) for each period. It then defines
a global support grid to ensure comparability across time and nterpolates density functions
onto this common grid. Finally, it creates a 3D surface plot (ridge plot) to visualize
distributional changes over time.

# Arguments

  - `data::Matrix{Float64}`: A matrix of density functions over time.

# Keyword Arguments

  - `n_grid_points::Int64`: Number of grid points for the common support grid (default
    `200`).
  - `legend::Bool`: If `true`, includes a legend in the plot (default `false`).

# Returns

  - A 3D surface plot of the IRF distributions over time.
"""
function plot_smooth_ridge(data::Matrix{Float64}; n_grid_points = 200, legend = false)

    # Compute KDE densities
    vec_data = [data[:, i] for i in axes(data, 2)]
    densities = [kde(dat) for dat in vec_data]
    global_min = minimum(minimum(den.x) for den in densities)
    global_max = maximum(maximum(den.x) for den in densities)

    # IRFs will change support over time -- for plotting, global grid needs to be defined
    common_grid = range(global_min, global_max; length = n_grid_points)

    # Interpolate densities onto common grid
    densities_interp = []
    for den in densities
        interp_func = interpolate((den.x,), den.density, Gridded(Linear()))
        extrapolated_func = extrapolate(interp_func, Line()) # Extrapolate with linear function
        push!(densities_interp, extrapolated_func.(common_grid)) # Interpolate onto common grid
    end

    # Place all densities togehter
    density_matrix = hcat(densities_interp...) # Now all KDEs are on the same x-grid
    x_axis = 1:size(density_matrix, 2) # Corresponds to the horizon (e.g., years)
    y_axis = sign.(common_grid) .* log.(abs.(common_grid) .+ 1)

    # Create 3D surface plot
    p = Plots.surface(
        x_axis,
        y_axis,
        density_matrix;
        camera = (60, 50),
        size = (600, 500),
        color = :winter,
        legend = legend,
        xlabel = "Horizon",
        ylabel = "Grid",
        zlabel = "Density",
        display_option = Plots.GR.OPTION_SHADED_MESH,
    )

    return p
end

"""
    plot_bivariate_animation(data; legend = false)

Creates a 3D animated visualization of how a bivariate distribution evolves over time.

# Arguments

  - `data::Array{Float64,3}`: A 3D array of bivariate distributions over time.

# Keyword Arguments

  - `legend::Bool`: If `true`, includes a legend in the plot (default `false`).

# Returns

  - An animated 3D surface plot of the bivariate distributions over time.
"""
function plot_bivariate_animation(data::Array{Float64,3}; legend = false)
    axes1 = 1:size(data, 1)
    axes2 = 1:size(data, 2)
    T = size(data, 3)

    anim = @animate for i ∈ 1:T
        Plots.surface(
            axes1,
            axes2,
            data[:, :, i];
            camera = (70, 30),
            size = (600, 500),
            color = :winter,
            legend = legend,
            display_option = Plots.GR.OPTION_SHADED_MESH,
        )
    end

    return anim
end
