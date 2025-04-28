"""
    compute_vardecomp(IRFs)

Computes the variance decomposition of the impulse response functions (IRFs) of a model.

# Arguments

  - `IRFs::Array{Float64}`: 3D array of IRFs for each exogenous variable, with dimensions,
    returned from `compute_irfs`:

      + 1: States and controls
      + 2: Time periods
      + 3: Exogenous variables

# Returns

  - `VDs`: 3D array of variance decompositions for each exogenous variable, with dimensions:

      + 1: States and controls
      + 2: Time periods
      + 3: Exogenous variables

  - `exovars_names`: Names of the exogenous variables.
"""
function compute_vardecomp(IRFs::Array{Float64})
    return cumsum(IRFs .^ 2.0; dims = 2) ./
           (sum(cumsum(IRFs .^ 2.0; dims = 2); dims = 3) .+ 10.0 * eps())
end

"""
    compute_vardecomp_bcfreq(exovars_full, stds, gx, hx; passband = (6, 32), ngrid = 512)

This function is designed to produce a variance decomposition at business cycle frequencies.
It produces a variance decomposition of the linearized solution. It returns the variance
decomposition at business cycle frequencies based on Uhlig (2001) and the unconditional
variance.

# Arguments

  - `exovars_full::Vector{Int64}`: Vector of positional indices of all exogenous variables
    to which shocks; note: here, no subset of exogenous variables is allowed.
  - `stds::Vector{Float64}`: Vector of standard deviations of the shocks.
  - `gx::Matrix{Float64}`: Control matrix (states to controls equations)
  - `hx::Matrix{Float64}`: Transition matrix for states (state transition equations)

# Keyword Arguments

  - `passband::Tuple{Int64,Int64}`: A tuple specifying the horizons associated with the
    business cycle. Default is (6, 32).
  - `ngrid::Int64`: The number of grid points for the computation of the band pass filter.
    Default is 512.

# Returns

  - `VD`: Variance decomposition at business cycle frequencies.
  - `UnconditionalVar`: Unconditional variance.
"""
function compute_vardecomp_bcfreq(
    exovars_full::Vector{Int64},
    stds::Vector{Float64},
    gx::Matrix{Float64},
    hx::Matrix{Float64};
    passband::Tuple{Int64,Int64} = (6, 32),
    ngrid::Int64 = 512,
)

    # Variables and indices
    nshocks = length(exovars_full)
    indexes_shocks = copy(exovars_full)
    nstates = size(hx, 1)
    indexes_states = collect(1:nstates)
    indexes_endostates = setdiff(indexes_states, indexes_shocks)
    nendostates = length(indexes_endostates)
    ncontrols = size(gx, 1)
    indexes_controls = collect((nstates + 1):(nstates + ncontrols))
    ntotal = ncontrols + nstates

    # Covariance matrix of states (and shocks)
    Σ = zeros(nstates, nstates)
    for (i, exovar) in enumerate(exovars_full)
        Σ[exovar, exovar] = stds[i]^2
    end
    Σ_exo = Σ[indexes_shocks, indexes_shocks]

    # Uhlig (2001) notation of matrices
    P = hx[indexes_endostates, indexes_endostates]
    Q = hx[indexes_endostates, indexes_shocks]
    R = gx[:, indexes_endostates]
    S = gx[:, indexes_shocks]
    N = hx[indexes_shocks, indexes_shocks]

    # Compute variance decomposition
    VD_alt, UnconditionalVar_alt =
        uhlig2001(P, Q, R, S, N, Σ_exo, ntotal, nendostates, nshocks, ngrid, passband)

    # Reordering in original format
    reordering = [indexes_endostates; indexes_controls; indexes_shocks]
    VD = similar(VD_alt)
    VD[reordering, :] = VD_alt
    UnconditionalVar = similar(UnconditionalVar_alt)
    UnconditionalVar[reordering] = UnconditionalVar_alt

    return VD, UnconditionalVar
end

function uhlig2001(P, Q, R, S, N, Σ_exo, ntotal, nendostates, nshocks, ngrid, passband)
    freqs = 0:((2 * pi) / ngrid):(2 * pi * (1 - 0.5 / ngrid))
    tpos = exp.(im .* freqs)
    tneg = exp.(-im .* freqs)

    filter_gain = zeros(1, ngrid)
    lowest_periodicity = passband[2]
    highest_periodicity = passband[1]
    highest_periodicity = max(2, highest_periodicity)

    # Set filter gain for frequencies within the passband
    passband_indices =
        (freqs .>= (2.0 * pi / lowest_periodicity)) .&
        (freqs .<= (2.0 * pi / highest_periodicity))
    filter_gain[passband_indices] .= 1.0

    # Set filter gain for negative frequencies within the passband
    negative_passband_indices =
        (freqs .<= (-2.0 * pi / lowest_periodicity + 2.0 * pi)) .&
        (freqs .>= (-2.0 * pi / highest_periodicity + 2.0 * pi))
    filter_gain[negative_passband_indices] .= 1.0

    # Make sure that shock covariance matrix is positive definite
    Σ_exo = Σ_exo + eps() * I

    # Containers for matrices
    VD = zeros(ntotal, nshocks)
    mat_bp_col = zeros(ComplexF64, nshocks + 1, ngrid, ntotal)
    mat1n = zeros(ComplexF64, nendostates + nshocks, nshocks)
    mat1p = zeros(ComplexF64, nshocks, nendostates + nshocks)
    mat2n = zeros(ComplexF64, nshocks, nshocks)
    mat2p = zeros(ComplexF64, nshocks, nshocks)
    mat3n = zeros(ComplexF64, ntotal, nendostates + nshocks)
    mat3p = zeros(ComplexF64, nendostates + nshocks, ntotal)
    mat31n = zeros(ComplexF64, ntotal, nshocks)
    mat13p = zeros(ComplexF64, nshocks, ntotal)
    inner = zeros(ComplexF64, nshocks, nshocks)
    f_omega = (1 / (2 * pi)) * inner
    g_omega = zeros(ComplexF64, ntotal, ntotal)
    f_bp = zeros(ComplexF64, ntotal, ntotal)

    for ig = 1:ngrid
        if filter_gain[ig] == 0.0
            mat_bp_col[:, ig, :] = zeros(ComplexF64, nshocks + 1, ntotal)
        else
            mat1n .= [(I - P * tneg[ig]) \ Q; I]
            mat1p .= [Q' / (I - P' * tpos[ig]) I]
            mat2n .= (I - N * tneg[ig])
            mat2p .= (I - N' * tpos[ig])
            mat3n .=
                [I zeros(nendostates, nshocks); R*tneg[ig] S; zeros(nshocks, nendostates) I]
            mat3p .= [
                I R'*tpos[ig] zeros(nendostates, nshocks)
                zeros(nshocks, nendostates) S' I
            ]
            mat31n .= mat3n * mat1n
            mat13p .= mat1p * mat3p
            inner .= ((mat2n \ Σ_exo) / mat2p)
            f_omega .= (1 / (2 * pi)) * inner
            g_omega .= mat31n * f_omega * mat13p
            f_bp .= filter_gain[ig]^2 * g_omega
            mat_bp_col[1, ig, :] = diag(f_bp)[:]
            for i = 1:nshocks
                Sigma_exo_i = zeros(size(Σ_exo))
                Sigma_exo_i[i, i] = Σ_exo[i, i]
                f_omega .= (1 / (2 * pi)) * ((mat2n \ Sigma_exo_i) / mat2p)
                g_omega .= mat31n * f_omega * mat13p
                f_bp .= filter_gain[ig]^2 * g_omega
                mat_bp_col[i + 1, ig, :] = diag(f_bp)[:]
            end
        end
    end

    imat_bp_col = real.(ifft(mat_bp_col, 2)) * (2 * pi)
    UnconditionalVar = imat_bp_col[1, 1, :][:]
    for i = 1:nshocks
        VD[:, i] = abs.(imat_bp_col[i + 1, 1, :]) ./ UnconditionalVar
    end

    return VD, UnconditionalVar
end
