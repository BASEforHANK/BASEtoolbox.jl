@doc raw"""
    Ksupply(RB_guess,R_guess,w_guess,profit_guess,n_par,m_par)

Calculate the aggregate savings when households face idiosyncratic income risk.

Idiosyncratic state is tuple ``(m,k,y)``, where
``m``: liquid assets, ``k``: illiquid assets, ``y``: labor income

# Arguments
- `R_guess`: real interest rate illiquid assets
- `RB_guess`: nominal rate on liquid assets
- `w_guess`: wages
- `profit_guess`: profits
- `n_par::NumericalParameters`
- `m_par::ModelParameters`

# Returns
- `K`,`B`: aggregate saving in illiquid (`K`) and liquid (`B`) assets
-  `TransitionMat`,`TransitionMat_a`,`TransitionMat_n`: `sparse` transition matrices
    (average, with [`a`] or without [`n`] adjustment of illiquid asset)
- `distr`: ergodic steady state of `TransitionMat`
- `c_a_star`,`m_a_star`,`k_a_star`,`c_n_star`,`m_n_star`: optimal policies for
    consumption [`c`], liquid [`m`] and illiquid [`k`] asset, with [`a`] or
    without [`n`] adjustment of illiquid asset
- `V_m`,`V_k`: marginal value functions
"""
function Ksupply(RB_guess::Float64, R_guess::Float64, n_par::NumericalParameters,
    m_par::ModelParameters, Vm::AbstractArray, Vk::AbstractArray, distr_guess::AbstractArray,
    inc::AbstractArray, eff_int::AbstractArray)

    #   initialize distance variables
    dist                = 9999.0
    dist1               = dist
    dist2               = dist

    q                   = 1.0       # price of Capital
    #----------------------------------------------------------------------------
    # Iterate over consumption policies
    #----------------------------------------------------------------------------
    count    = 0
    n        = size(Vm)
    # containers for policies, marginal value functions etc.
    m_n_star = similar(Vm);         m_a_star = similar(Vm)
    k_a_star = similar(Vm)
    c_a_star = similar(Vm);         c_n_star = similar(Vm)
    EVm      = similar(Vm);         EVk      = similar(Vk)
    Vm_new   = similar(Vm);         Vk_new   = similar(Vk)    
    iVm      = invmutil(Vm,m_par);  iVk      = invmutil(Vk,m_par)
    iVm_new  = similar(iVm);        iVk_new  = similar(iVk)
    E_return_diff = similar(EVm);   EMU      = similar(EVm)
    c_star_n      = similar(EVm);   m_star_n = similar(EVm)
    mutil_c_a     = similar(EVm);   D1       = similar(EVm)
    D2       = similar(EVm)
    Resource_grid = reshape(inc[2].+ inc[3] .+ inc[4] ,(n[1].*n[2], n[3]))
    while dist > n_par.ϵ && count < 10000 # Iterate consumption policies until converegence
        count           = count + 1
        # Take expectations for labor income change
        #EVk  .= reshape(reshape(Vk, (n[1] .* n[2], n[3])) * n_par.Π', (n[1], n[2], n[3]))
        BLAS.gemm!('N', 'T', 1.0, reshape(Vk, (n[1] .* n[2], n[3])), n_par.Π, 0.0, reshape(EVk, (n[1] .* n[2], n[3])))
        EVk .= reshape(EVk,(n[1], n[2], n[3]))
        BLAS.gemm!('N', 'T', 1.0, reshape(Vm, (n[1] .* n[2], n[3])), n_par.Π, 0.0, reshape(EVm, (n[1] .* n[2], n[3])))
        EVm .= reshape(EVm,(n[1], n[2], n[3]))
        EVm .*= eff_int

        # Policy update step
        EGM_policyupdate!(c_a_star, m_a_star, k_a_star, c_n_star, m_n_star,
                          E_return_diff, EMU, c_star_n, m_star_n, Resource_grid, 
                          EVm, EVk, q, m_par.π, RB_guess, 1.0, inc, n_par, m_par, false)

        # marginal value update step
        updateV!(Vk_new, Vm_new, mutil_c_a, EVk, c_a_star, c_n_star, m_n_star, R_guess - 1.0, q, m_par, n_par)
        invmutil!(iVk_new, Vk_new, m_par)
        invmutil!(iVm_new, Vm_new, m_par)
        # Calculate distance in updates
        D1            .= iVk_new .- iVk
        D2            .= iVm_new .- iVm
        dist1          = maximum(abs, D1)
        dist2          = maximum(abs, D2)
        dist           = max(dist1, dist2) # distance of old and new policy

        # update policy guess/marginal values of liquid/illiquid assets
        Vm             .= Vm_new
        Vk             .= Vk_new
        iVk            .= iVk_new
        iVm            .= iVm_new
    end
    println("EGM Iterations: ", count)    
    println("EGM Dist: ", dist)    

    #------------------------------------------------------
    # Find stationary distribution (Is direct transition better for large model?)
    #------------------------------------------------------

    # Define transition matrix
    S_a, T_a, W_a, S_n, T_n, W_n    = MakeTransition(m_a_star,  m_n_star, k_a_star, n_par.Π, n_par)
    TransitionMat_a                 = sparse(S_a, T_a, W_a, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    TransitionMat_n                 = sparse(S_n, T_n, W_n, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    Γ                               = m_par.λ .* TransitionMat_a .+ (1.0 .- m_par.λ) .* TransitionMat_n

    # Calculate left-hand unit eigenvector
    
    aux = real.(eigsolve(Γ', distr_guess[:], 1)[2][1])
   
    ## Exploit that the Eigenvector of eigenvalue 1 is the nullspace of TransitionMat' -I
    #     Q_T = LinearMap((dmu, mu) -> dist_change!(dmu, mu, Γ), n_par.nm * n_par.nk * n_par.ny, ismutating = true)
    #     aux = fill(1/(n_par.nm * n_par.nk * n_par.ny), n_par.nm * n_par.nk * n_par.ny)#distr_guess[:] # can't use 0 as initial guess
    #     gmres!(aux, Q_T, zeros(n_par.nm * n_par.nk * n_par.ny))  # i.e., solve x'(Γ-I) = 0 iteratively
    ##qr algorithm for nullspace finding
    #     aux2 = qr(Γ - I)
    #     aux = Array{Float64}(undef, n_par.nm * n_par.nk * n_par.ny)
    #     aux[aux2.prow] = aux2.Q[:,end]
    #
    distr = (reshape((aux[:]) ./ sum((aux[:])),  (n_par.nm, n_par.nk, n_par.ny)))

    #-----------------------------------------------------------------------------
    # Calculate capital stock
    #-----------------------------------------------------------------------------
    K = dot(sum(distr,dims=(1,3)), n_par.grid_k)
    B = dot(sum(distr,dims=(2,3)), n_par.grid_m)
    return K, B, Γ, TransitionMat_a, TransitionMat_n, c_a_star, m_a_star, k_a_star, c_n_star, m_n_star, Vm, Vk, distr
end
function next_dist(mu,Q)
    @unpack m,n, colptr, rowval, nzval = Q
    nextmu=similar(mu)
    @inbounds for col = 1:m
        nextmu[col] = 0.0
        for n_row = colptr[col]:colptr[col+1]-1
            nextmu[col] += mu[rowval[n_row]]*nzval[n_row]
        end
    end
    return nextmu
end