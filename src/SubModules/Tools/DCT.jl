"""
    mydctmx(n::Int)

Generate the Discrete Cosine Transform (DCT) matrix of size `n × n`.

This implementation uses the DCT-II, i.e., it normalizes the matrix to be orthonormal.

# Arguments

  - `n::Int`: The size of the square DCT matrix, corresponding to the number of grid points
    of the to be transformed dimension.

# Returns

  - `Array{Float64, 2}`: A square `n × n` matrix representing the Discrete Cosine Transform.
"""
function mydctmx(n::Int)
    DC::Array{Float64,2} = zeros(n, n)
    for j = 0:(n - 1)
        DC[1, j + 1] = float(1 / sqrt(n))
        for i = 1:(n - 1)
            DC[i + 1, j + 1] = (pi * float((j + 1 / 2) * i / n))
            DC[i + 1, j + 1] = sqrt(2 / n) .* cos.(DC[i + 1, j + 1])
        end
    end
    return DC
end

"""
    uncompress(compressionIndexes, XC, DC, IDC)

Reconstructs the full 3D array from its compressed representation `XC`. Hence its the
inverse operation to [`Tools.compress()`](@ref).

The function first reverses the compression by replacing non-retained coefficients with
zeros, where the indices of the retained coefficients are mapped back to the original grid
via `compressionIndexes`. It then reverts the 3D DCT using the transformation matrices `DC`
and their corresponding inverse `IDC`.

# Arguments

  - `compressionIndexes::Vector{Int}`: Indices mapping the compressed data to the
    uncompressed data so that index `j` in `compressionIndexes` maps to the `j`-th most
    important index in the uncompressed data.

  - `XC::AbstractArray`: Compressed data as a 1D array.
  - `DC::AbstractArray`: Array of DCT matrices for each dimension:

      + `DC[1]`: DCT matrix for the first dimension (illiquid asset grid).
      + `DC[2]`: DCT matrix for the second dimension (liquid asset grid).
      + `DC[3]`: DCT matrix for the third dimension (productivity grid).
  - `IDC::AbstractArray`: Array of inverse DCT matrices for each dimension:

      + `IDC[1]`: Inverse DCT matrix for the first dimension (illiquid asset grid).
      + `IDC[2]`: Inverse DCT matrix for the second dimension (liquid asset grid).
      + `IDC[3]`: Inverse DCT matrix for the third dimension (productivity grid).

# Returns

  - `θ::Vector`: Reconstructed (uncompressed) data, reshaped into a 1D array, i.e. all 3
    dimensions are stacked.

# Note

  - Since usually information is lost during compression, the uncompressed data will not be
    identical to the original data.
"""
function uncompress(compressionIndexes, XC, DC, IDC)
    nb = size(DC[1], 1)
    nk = size(DC[2], 1)
    nh = size(DC[3], 1)
    θ1 = zeros(eltype(XC), nb, nk, nh)
    for j in eachindex(XC)
        θ1[compressionIndexes[j]] = copy(XC[j])
    end
    @views for hh = 1:nh
        θ1[:, :, hh] = IDC[1] * θ1[:, :, hh] * DC[2]
    end
    @views for bb = 1:nb
        θ1[bb, :, :] = θ1[bb, :, :] * DC[3]
    end
    θ = reshape(θ1, (nb) * (nk) * (nh))
    return θ
end

"""
    compress(compressionIndexes, XU, DC, IDC)

Compresses the array `XU` using a 3D Discrete Cosine Transform (DCT) and retains only the
most relevant coefficients.

The function first applies the 3D DCT using the transformation matrices `DC` and their
corresponding inverse `IDC`. It then selects a subset of the transformed coefficients based
on `compressionIndexes`.

# Arguments

  - `compressionIndexes::AbstractArray`: Indices specifying which coefficients of the
    transformed data should be retained.

  - `XU::AbstractArray`: The input array to be compressed, `(nb, nk, nh)`.
  - `DC::AbstractArray`: Array of DCT matrices for each dimension:

      + `DC[1]`: DCT matrix for the first dimension (illiquid asset grid).
      + `DC[2]`: DCT matrix for the second dimension (liquid asset grid).
      + `DC[3]`: DCT matrix for the third dimension (productivity grid).
  - `IDC::AbstractArray`: Array of inverse DCT matrices for each dimension:

      + `IDC[1]`: Inverse DCT matrix for the first dimension (illiquid asset grid).
      + `IDC[2]`: Inverse DCT matrix for the second dimension (liquid asset grid).
      + `IDC[3]`: Inverse DCT matrix for the third dimension (productivity grid).

# Returns

  - `θ::AbstractArray`: Compressed data as a 1D array containing only the selected
    coefficients.
"""
function compress(
    compressionIndexes::AbstractArray,
    XU::AbstractArray,
    DC::AbstractArray,
    IDC::AbstractArray,
)
    nb, nk, nh = size(XU)
    θ = zeros(eltype(XU), length(compressionIndexes))
    XU2 = similar(XU)
    @inbounds @views for b = 1:nb
        XU2[b, :, :] = DC[2] * XU[b, :, :] * IDC[3]
    end
    @inbounds @views for y = 1:nh
        XU2[:, :, y] = DC[1] * XU2[:, :, y]
    end
    θ = XU2[compressionIndexes]
    return θ
end

"""
    select_ind(Theta, reduc_value::Float64)

Select the multidimensional DCT coefficients that maintain a specified fraction of the
"energy" of `Theta`, (`1 - reduc_value`) of the total energy.

The function iteratively selects coefficients from the Discrete Cosine Transformation (DCT)
`Theta`, ensuring that the retained coefficients preserve a certain percentage of the total
energy of `Theta`. The coefficients are sorted by their absolute size, and the function
stops when the retained coefficients maintain at least `1 - reduc_value` of the energy.

This share of "energy" is equivalent to the share of the Euclidean norm of `Theta` explained
by the retained coefficients.

# Arguments

  - `Theta::AbstractArray`: The DCT-transformed array, typically representing marginal
    liquid asset values or more generally, a value function.
  - `reduc_value::Float64`: The fraction of the energy to be discarded from the original
    `Theta`. The retained coefficients preserve `1 - reduc_value` of the total energy.

# Returns

  - `select_ind::AbstractArray`: Indices of the retained coefficients that preserve the
    desired energy fraction.
"""
function select_ind(Theta, reduc_value::Float64)
    Theta = Theta[:]

    # Indexes of coefficients sorted by their absolute size
    ind = sortperm(abs.(Theta[:]); rev = true)

    # Container to store the number of retained coefficients
    coeffs = 1

    # Find the important basis functions (discrete cosine) for Theta and add retained until
    # only n_par.reduc_value share of energy is lost
    while norm(Theta[ind[1:coeffs]]) / norm(Theta) < 1 - reduc_value
        coeffs += 1
    end

    # The indices of the retained coefficients
    select_ind = ind[1:coeffs]

    return select_ind
end
