##########################################################
# Matrix to remove one degree of freedom from distribution
#---------------------------------------------------------
function shuffleMatrix(distr, nb, nk, nh)
    distr_b = sum(sum(distr; dims = 3); dims = 2) ./ sum(distr[:])
    distr_k = sum(sum(distr; dims = 3); dims = 1) ./ sum(distr[:])
    distr_h = sum(sum(distr; dims = 2); dims = 1) ./ sum(distr[:])
    Γ = Array{Array{Float64,2},1}(undef, 3)
    Γ[1] = zeros(Float64, (nb, nb - 1))
    Γ[2] = zeros(Float64, (nk, nk - 1))
    Γ[3] = zeros(Float64, (nh, nh - 1))
    for j = 1:(nb - 1)
        Γ[1][:, j] = -distr_b[:]
        Γ[1][j, j] = 1 - distr_b[j]
        Γ[1][j, j] = Γ[1][j, j] - sum(Γ[1][:, j])
    end
    for j = 1:(nk - 1)
        Γ[2][:, j] = -distr_k[:]
        Γ[2][j, j] = 1 - distr_k[j]
        Γ[2][j, j] = Γ[2][j, j] - sum(Γ[2][:, j])
    end
    for j = 1:(nh - 1)
        Γ[3][:, j] = -distr_h[:]
        Γ[3][j, j] = 1 - distr_h[j]
        Γ[3][j, j] = Γ[3][j, j] - sum(Γ[3][:, j])
    end

    return Γ
end
