
@doc raw"""
    nearest_spd(A)

Return the nearest (in Frobenius norm) Symmetric Positive Definite matrix to `A`.

Based on [answer in MATLAB Central forum](https://de.mathworks.com/matlabcentral/answers/320134-make-sample-covariance-correlation-matrix-positive-definite).
From Higham: "The nearest symmetric positive semidefinite matrix in the
Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
where H is the symmetric polar factor of B=(A + A')/2."

# Arguments
`A`: square matrix

# Returns
`Ahat`: nearest SPD matrix to `A`
"""
function nearest_spd(A)

    # symmetrize A into B
    B = 0.5 .* (A .+ A')
    FU, FS, FVt = LinearAlgebra.LAPACK.gesvd!('N', 'S', copy(B))
    H = FVt' * Diagonal(FS) * FVt

    # get Ahat in the above formula
    Ahat = 0.5 .* (B .+ H)

    # ensure symmetry
    Ahat .= 0.5 .* (Ahat .+ Ahat')

    # test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    p = false
    k = 0
    count = 1
    while p == false && count < 100
        R = cholesky(Ahat; check = false)
        k += 1
        count = count + 1
        if ~issuccess(R)
            # Ahat failed the chol test. It must have been just a hair off,
            # due to floating point trash, so it is simplest now just to
            # tweak by adding a tiny multiple of an identity matrix.
            mineig = eigmin(Ahat)
            Ahat .+= (-mineig .* k .^ 2 .+ eps(mineig)) .* Diagonal(ones(size(Ahat, 1)))
        else
            p = true
        end
    end

    return Ahat
end
