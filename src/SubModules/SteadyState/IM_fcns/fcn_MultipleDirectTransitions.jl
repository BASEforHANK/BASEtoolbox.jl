@doc raw"""
    ...
"""
function MultipleDirectTransition(
    b_a_star::Array{Float64,3},
    b_n_star::Array{Float64,3},
    k_a_star::Array{Float64,3},
    distr::Array{Float64,3},
    λ::Float64,
    Π::Array{Float64,2},
    n_par,
)
    idk_a, wR_k_a, wL_k_a = MakeWeights(k_a_star, n_par.grid_k)
    idb_a, wR_b_a, wL_b_a = MakeWeights(b_a_star, n_par.grid_b)
    idb_n, wR_b_n, wL_b_n = MakeWeights(b_n_star, n_par.grid_b)
    dist = 1.0
    count = 1
    blockindex = (0:(n_par.nh - 1)) * n_par.nk * n_par.nb
    while (dist > n_par.ϵ) && (count < 10000)
        dPrime = zeros(typeof(distr[1]), size(distr))
        for hh = 1:(n_par.nh) # all current income states
            for kk = 1:(n_par.nk) # all current illiquid asset states
                #idk_n = kk
                for bb = 1:(n_par.nb)
                    dd = distr[bb, kk, hh]
                    IDD_a = idb_a[bb, kk, hh] .+ (idk_a[bb, kk, hh] - 1) .* n_par.nb
                    IDD_n = idb_n[bb, kk, hh] .+ (kk - 1) .* n_par.nb
                    DLL_a = dd .* wL_k_a[bb, kk, hh] .* wL_b_a[bb, kk, hh]
                    DLR_a = dd .* wL_k_a[bb, kk, hh] .* wR_b_a[bb, kk, hh]
                    DRL_a = dd .* wR_k_a[bb, kk, hh] .* wL_b_a[bb, kk, hh]
                    DRR_a = dd .* wR_k_a[bb, kk, hh] .* wR_b_a[bb, kk, hh]
                    DL_n = dd .* wL_b_n[bb, kk, hh]
                    DR_n = dd .* wR_b_n[bb, kk, hh]
                    for hh_prime = 1:(n_par.nh)
                        id_a = IDD_a + blockindex[hh_prime]
                        id_n = IDD_n + blockindex[hh_prime]
                        fac = λ .* Π[hh, hh_prime]
                        dPrime[id_a] += fac .* DLL_a
                        dPrime[id_a + 1] += fac .* DLR_a
                        dPrime[id_a + n_par.nb] += fac .* DRL_a
                        dPrime[id_a + n_par.nb + 1] += fac .* DRR_a
                        dPrime[id_n] += (1.0 - λ) .* Π[hh, hh_prime] .* DL_n
                        dPrime[id_n + 1] += (1.0 - λ) .* Π[hh, hh_prime] .* DR_n
                    end
                end
            end
        end
        dist = maximum(abs.(dPrime[:] - distr[:]))
        distr = dPrime
        count = count + 1
    end
    return distr, dist, count
end
