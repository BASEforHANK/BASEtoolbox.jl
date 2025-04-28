"""
    compute_irfs(
      exovars,
      gx,
      hx,
      XSS,
      ids;
      T = 1000,
      init_val = fill(0.01, length(exovars)),
      verbose = true
      distribution = false
      comp_ids = nothing
    )

Computes impulse response functions (IRFs) for a given set of shocks to exogenous variables.

# Arguments

  - `exovars::Vector{Int64}`: Vector of positional indices of exogenous variables to which
    shocks are applied.
  - `gx::Matrix{Float64}`: Control matrix (states to controls equations)
  - `hx::Matrix{Float64}`: Transition matrix for states (state transition equations)
  - `XSS::Vector{Float64}`: Steady state values of the model variables.
  - `ids`: Indexes of the model variables.

# Keyword Arguments

  - `T::Int64`: Number of periods for which to compute IRFs.
  - `init_val::Vector{Float64}`: Initial value of the shock to each exogenous variable,
    defaults to 0.01 for all shocks.
  - `verbose::Bool`: Print progress to console.
  - `distribution::Bool`: Compute distributional IRFs, defaults to false.
  - `comp_ids`: Compression indices for the distribution, as created by
    `prepare_linearization`. Needed if `distribution` is true. Defaults to nothing.

# Returns

  - `IRFs`, `IRFs_lvl`: 3D array of (level) IRFs for each exogenous variable, with
    dimensions:

      + 1: States and controls
      + 2: Time periods
      + 3: Exogenous variables

  - `IRFs_order`: Names of the exogenous variables and their indices in the IRFs array.
  - `IRFs_dist`: Dictionary containing the distributional IRFs in levels, if requested.
    Dimensions are the same as for `IRFs` and `IRFs_lvl`, with the following keys:

      + "Wb": Marginal value function of liquid assets, nb x T.
      + "Wk": Marginal value function of illiquid assets, nk x T.
      + "distr_b": Distribution of liquid assets, nb x T.
      + "distr_k": Distribution of illiquid assets, nk x T.
      + "distr": Joint distribution of liquid and illiquid assets, nb x nk x T.
"""
function compute_irfs(
    exovars::Vector{Int64},
    gx::Matrix{Float64},
    hx::Matrix{Float64},
    XSS::Vector{Float64},
    ids;
    T::Int64 = 1000,
    init_val::Vector{Float64} = fill(0.01, length(exovars)),
    verbose::Bool = true,
    distribution::Bool = false,
    comp_ids = nothing,
)

    # If distributional IRFs are requested, check if the compression indices are provided
    if distribution && isnothing(comp_ids)
        throw(ArgumentError("Missing `comp_ids`: Required for distributional IRFs."))
    end

    # Compute the number of states and controls from gx and hx
    ncontrols = size(gx, 1)
    nstates = size(hx, 1)

    # Initialize IRFs for all selected exogenous variables
    IRFs = zeros(nstates + ncontrols, T, length(exovars))
    IRFs_lvl = similar(IRFs)
    if distribution
        IRFs_dist = Dict(
            "Wb" => zeros(length(ids.distr_bSS), T, length(exovars)),
            "Wk" => zeros(length(ids.distr_kSS), T, length(exovars)),
            "distr_b" => zeros(length(ids.distr_bSS), T, length(exovars)),
            "distr_k" => zeros(length(ids.distr_kSS), T, length(exovars)),
            "distr" =>
                zeros(length(ids.distr_bSS), length(ids.distr_kSS), T, length(exovars)),
        )
    end

    # Store the shock names
    IRFs_order = [find_field_with_value(ids, exovars[i], false) for i = 1:length(exovars)]

    # Compute IRFs for each exogenous variable
    for (i, exovar) in enumerate(exovars)
        if verbose
            @printf "Computing IRFs for %s with index %d and initial condition of %f.\n" IRFs_order[i] exovar init_val[i]
        end
        if distribution
            IRFs[:, :, i],
            IRFs_lvl[:, :, i],
            IRFs_dist["Wb"][:, :, i],
            IRFs_dist["Wk"][:, :, i],
            IRFs_dist["distr_b"][:, :, i],
            IRFs_dist["distr_k"][:, :, i],
            IRFs_dist["distr"][:, :, :, i] = compute_irfs_inner(
                exovars[i],
                gx,
                hx,
                XSS,
                ids,
                T,
                nstates,
                ncontrols,
                init_val[i],
                distribution,
                comp_ids,
            )
        else
            IRFs[:, :, i], IRFs_lvl[:, :, i] = compute_irfs_inner(
                exovars[i],
                gx,
                hx,
                XSS,
                ids,
                T,
                nstates,
                ncontrols,
                init_val[i],
                distribution,
                comp_ids,
            )
        end
    end

    if distribution
        return IRFs, IRFs_lvl, IRFs_order, IRFs_dist
    else
        return IRFs, IRFs_lvl, IRFs_order
    end
end

"""
    compute_irfs_inner(
        exovar,
        gx,
        hx,
        XSS,
        ids,
        T,
        nstates,
        ncontrols,
        init_val,
        distribution,
        comp_ids
    )

Computes impulse response functions (IRFs) for a given shock to a single exogenous variable.

See `compute_irfs` for argument and return value descriptions.
"""
function compute_irfs_inner(
    exovar::Int64,
    gx::Matrix{Float64},
    hx::Matrix{Float64},
    XSS::Vector{Float64},
    ids,
    T::Int64,
    nstates::Int64,
    ncontrols::Int64,
    init_val::Float64,
    distribution::Bool,
    comp_ids,
)

    # Initialize matrices for states and controls
    S_t = zeros(nstates, T)
    C_t = zeros(ncontrols, T)

    # Initial conditions: states by assumption, controls as implied by gx and initial state
    S_t[exovar, 1] = init_val
    C_t[:, 1] = gx * S_t[:, 1]

    # Simulation: iterate forward
    for t = 2:T
        S_t[:, t] = hx * S_t[:, t - 1]
        C_t[:, t] = gx * S_t[:, t]
    end

    # Recompute levels for the original IRFs, as defined in macro @generate_equations
    original = [S_t; C_t]
    level = fill(NaN64, size(original))

    # Start with the aggregate variables
    idx = [getfield(ids, Symbol(j)) for j in aggr_names]
    idxSS = [getfield(ids, Symbol(j, "SS")) for j in aggr_names]
    level[idx, :] = exp.(XSS[idxSS] .+ original[idx, :])

    if distribution
        WbIRF, WkIRFs, distr_bIRF, distr_kIRF, distrIRF =
            compute_irfs_inner_distribution(original, ids, XSS, comp_ids)
        return original, level, WbIRF, WkIRFs, distr_bIRF, distr_kIRF, distrIRF
    else
        return original, level
    end
end

"""
    compute_irfs_inner_distribution(original, ids, XSS, comp_ids)

Compute impulse response functions (IRFs) for the distribution and marginal value functions
for a given shock to a single exogenous variable.

This function reconstructs the impulse responses of wealth distributions and marginal value
functions from the compressed IRFs. It first extracts steady-state values and applies the
discrete cosine transform (DCT) to uncompress the shocks. It then perturbs the steady-state
distributions with the IRFs and computes the updated marginal value functions and joint
distribution of liquid and illiquid assets.

The function uses a precomputed shuffle matrix to adjust distributions appropriately before
summing over dimensions to extract relevant marginal values.

# Arguments

  - `original::Array{Float64,2}`: The original IRF data matrix containing the compressed
    IRFs for the given exogenous variable.
  - `ids`: Indexes of the model variables.
  - `XSS::Vector{Float64}`: Steady state values of the model variables.
  - `comp_ids`: Compression indices for the distribution.

# Returns

  - `WbIRF::Matrix{Float64}`: The impulse responses for the marginal value function of
    liquid assets over time, size `nb x T`.
  - `WkIRFs::Matrix{Float64}`: The impulse responses for the marginal value function of
    illiquid assets over time, size `nk x T`.
  - `distr_bIRF::Matrix{Float64}`: The impulse responses for the distribution of liquid
    assets over time, size `nb x T`.
  - `distr_kIRF::Matrix{Float64}`: The impulse responses for the distribution of illiquid
    assets over time, size `nk x T`.
  - `distrIRF::Matrix{Float64}`: The impulse responses for the joint distribution of liquid
    and illiquid assets, size `nb x nk x T`.
"""
function compute_irfs_inner_distribution(
    original::Array{Float64,2},
    ids,
    XSS::Vector{Float64},
    comp_ids,
)

    # Grid sizes
    nb = length(ids.distr_bSS)
    nk = length(ids.distr_kSS)
    nh = length(ids.distr_hSS)
    T = size(original, 2)

    # Get steady-state objects
    WbSS = XSS[ids.WbSS]
    WkSS = XSS[ids.WkSS]
    COPSS = XSS[ids.COPSS]
    distr_bSS = XSS[ids.distr_bSS]
    distr_kSS = XSS[ids.distr_kSS]
    distr_SS = reshape(XSS[ids.COPSS], (nb, nk, nh))

    # DCT containers
    DC = Array{Array{Float64,2},1}(undef, 3)
    DC[1] = mydctmx(nb)
    DC[2] = mydctmx(nk)
    DC[3] = mydctmx(nh)
    IDC = [DC[1]', DC[2]', DC[3]']

    # Shuffle matrix
    Γ = shuffleMatrix(distr_SS, nb, nk, nh)

    # Define IRF containers
    WbIRF = zeros(nb, T)
    WkIRFs = zeros(nk, T)
    distr_bIRF = zeros(nb, T)
    distr_kIRF = zeros(nk, T)
    distrIRF = zeros(nb, nk, T)

    # Compute IRFs
    for t = 1:T

        # Uncompress IRFs + perturb steady-state values with IRFs
        Wb_full = exp.(WbSS .+ uncompress(comp_ids[1], original[ids.Wb, t], DC, IDC))
        Wk_full = exp.(WkSS .+ uncompress(comp_ids[2], original[ids.Wk, t], DC, IDC))
        CP_full = exp.(COPSS .+ uncompress(comp_ids[3], original[ids.COP, t], DC, IDC))

        # Masses of liquid and illiquid assets
        distr_bIRF[:, t] = distr_bSS .+ Γ[1] * original[ids.distr_b, t]
        distr_kIRF[:, t] = distr_kSS .+ Γ[2] * original[ids.distr_k, t]

        # Get the marginal value function of liquid assets
        WbIRF[:, t] = sum(reshape(Wb_full, (nb, nk, nh)); dims = (2, 3))[:]

        # Get the marginal value function of illiquid assets
        WkIRFs[:, t] = sum(reshape(Wk_full, (nb, nk, nh)); dims = (1, 3))[:]

        # Get the joint distribution of liquid and illiquid assets
        distrIRF[:, :, t] = sum(reshape(CP_full, (nb, nk, nh)); dims = 3)[:, :]
    end

    return WbIRF, WkIRFs, distr_bIRF, distr_kIRF, distrIRF
end
