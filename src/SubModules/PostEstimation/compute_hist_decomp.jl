"""
    compute_hist_decomp(exovars_full, gx, hx, smoother_output, ids)

Computes the historical decomposition of the model's variables in response to a set of
shocks, including the effect of initial conditions and each shock.

# Arguments

  - `exovars_full::Vector{Int64}`: Vector of positional indices of all exogenous variables
    to which shocks; note: here, no subset of exogenous variables is allowed.
  - `gx::Matrix{Float64}`: Control matrix (states to controls equations)
  - `hx::Matrix{Float64}`: Transition matrix for states (state transition equations)
  - `smoother_output`: The output from a smoother function.
  - `ids`: Indexes of the model variables.

# Returns

  - `ShockContr::Array{Float64}`: A 3D array representing the decomposition of the effects
    of each shock and the initial condition on states and controls over time. The dimensions
    are:

      + 1: States and controls
      + 2: Time periods
      + 3: Shocks (including an additional index for the initial condition)

  - `ShockContr_order::Vector{Any}`: A vector containing the names or identifiers of the
    shocks and the initial condition, in the order they appear in the `ShockContr` array.
"""
function compute_hist_decomp(
    exovars_full::Vector{Int64},
    gx::Matrix{Float64},
    hx::Matrix{Float64},
    smoother_output,
    ids,
)

    # Variables and indices
    nshocks = length(exovars_full)
    nstates = size(hx, 1)
    ncontrols = size(gx, 1)
    ntotal = ncontrols + nstates
    T = size(smoother_output[6], 2)

    ShockContr = Array{Float64}(undef, ntotal, T, nshocks + 1)

    # Compute effect of initial condition
    ShockContr_aux = zeros(ntotal, T)
    x = smoother_output[3][:, 1]
    for t = 1:T
        ShockContr_aux[:, t] = [I; gx] * x
        x = hx * x
    end
    ShockContr[:, :, nshocks + 1] = ShockContr_aux

    # Compute effect of each shock
    for j = 1:nshocks
        ShockContr_aux = zeros(ntotal, T)
        x = zeros(nstates)
        shock_index = exovars_full[j]
        for t = 1:T
            ShockContr_aux[:, t] = [I; gx] * x
            x = hx * x
            x[shock_index] += smoother_output[6][shock_index, t]
        end
        ShockContr[:, :, j] = ShockContr_aux
    end

    # Store order of shocks
    ShockContr_order =
        [find_field_with_value(ids, exovars_full[i], false) for i = 1:nshocks]
    ShockContr_order = [ShockContr_order; :initial]

    return ShockContr, ShockContr_order
end
