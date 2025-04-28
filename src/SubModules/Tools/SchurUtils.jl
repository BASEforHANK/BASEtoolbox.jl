"""
    real_schur(A, B)

Computes the real generalized Schur decomposition of the matrix pair `(A, B)`. Selects the
eigenvalues related to the state variables.

# Arguments

  - `A`: A square matrix.
  - `B`: A square matrix of the same size as `A`.

# Returns

  - `F`: The Schur decomposition of `(A, B)`, return of `LinearAlgebra.schur`.
  - `select_ev`: A boolean vector indicating which eigenvalues have a modulus greater than
    or equal to 1, i.e. which eigenvalues and/or eigenvectors are related to the state
    variables.
  - `nk`: The number of eigenvalues with a modulus greater than or equal to 1, i.e. the
    number of state variables.
  - `λ`: The vector of all generalized eigenvalues.
"""
function real_schur(A, B)
    F = LinearAlgebra.schur(A, B)
    α::Vector{complex(promote_type(eltype(A), eltype(B)))} = F.alpha
    λ = abs.(α) ./ abs.(F.beta)
    select_ev = λ .>= 1.0
    nk = sum(select_ev)
    return F, select_ev, nk, λ
end

"""
    complex_schur(A, B)

Version of [`BASEforHANK.Tools.real_schur()`](@ref) robust to complex input arguments. See
[`BASEforHANK.Tools.real_schur()`](@ref) for more details.
"""
function complex_schur(A, B)
    F = LinearAlgebra.schur(complex(A), complex(B))
    α::Vector{complex(promote_type(eltype(A), eltype(B)))} = F.alpha
    λ = abs.(α) ./ abs.(F.beta)
    select_ev = λ .>= 1.0
    nk = sum(select_ev)
    return F, select_ev, nk, λ
end
