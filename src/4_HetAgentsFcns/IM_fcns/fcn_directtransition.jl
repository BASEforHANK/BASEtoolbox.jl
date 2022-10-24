function DirectTransition(m_a_star::Array,
    m_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par::NumericalParameters)

    dPrime = zeros(eltype(distr),size(distr))
    DirectTransition!(dPrime,m_a_star, m_n_star,  k_a_star,  distr,  λ, Π, n_par)
    return dPrime
end
function DirectTransition!(dPrime,
    m_a_star::Array,
    m_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par::NumericalParameters)
    
    idk_a, wR_k_a = MakeWeightsLight(k_a_star,n_par.grid_k)
    idm_a, wR_m_a = MakeWeightsLight(m_a_star,n_par.grid_m)
    idm_n, wR_m_n = MakeWeightsLight(m_n_star,n_par.grid_m)
    blockindex = (0:n_par.ny-1)*n_par.nk*n_par.nm
    @inbounds begin
    for zz = 1:n_par.ny # all current income states
        for kk = 1:n_par.nk # all current illiquid asset states
            #idk_n = kk
            for mm = 1:n_par.nm
                dd    = distr[mm,kk,zz]
                IDD_a = (idm_a[mm,kk,zz].+(idk_a[mm,kk,zz] .-1).*n_par.nm)
                IDD_n = (idm_n[mm,kk,zz].+(kk-1).*n_par.nm)
                # liquid assets of non adjusters
                w     = wR_m_n[mm,kk,zz]
                DL_n  = (1.0 .- λ).*(dd.*(1.0 .- w))
                DR_n  = (1.0 .- λ).*(dd.*w)
                # illiquid assets of adjusters
                w     = wR_k_a[mm,kk,zz]
                dl    = λ.*(dd.*(1.0 .- w ))
                dr    = λ.*(dd.*w)
                # liquid assets of adjusters
                w     = wR_m_a[mm,kk,zz]
                DLL_a = (dl.*(1.0 .- w))
                DLR_a = (dl.*w)
                DRL_a = (dr.*(1.0 .- w))
                DRR_a = (dr.*w)
                for yy = 1:n_par.ny # add income transitions
                    pp   = Π[zz,yy]
                    id_a = IDD_a .+ blockindex[yy]
                    id_n = IDD_n .+ blockindex[yy]
                    dPrime[id_a]            += pp.*DLL_a
                    dPrime[id_a+1]          += pp.*DLR_a
                    dPrime[id_a+n_par.nm]   += pp.*DRL_a
                    dPrime[id_a+n_par.nm+1] += pp.*DRR_a
                    dPrime[id_n]            += pp.*DL_n
                    dPrime[id_n+1]          += pp.*DR_n
                end
            end
        end
    end
end

end

