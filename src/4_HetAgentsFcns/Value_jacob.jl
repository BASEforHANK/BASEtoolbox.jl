function Value_jacob_dct(VkSS, VmSS, n_par, m_par, mcw, A, q, RL, τprog, τlev, H, Ht, 
                            π, r, w, N, profits, unionprofits, av_tax_rate, σ)

    price = [mcw A q RL τprog τlev H Ht π r w N profits unionprofits av_tax_rate σ ][:]
    Γ,Vk,Vm = Transition(price, VkSS, VmSS, n_par, m_par)

    ff(p) = Values(p, VkSS, VmSS,Γ, n_par, m_par)

    J    = ForwardDiff.jacobian(p->ff(p), price)
    println("Jacobian done!")
    N    = n_par.nm*n_par.nk*n_par.ny
    Jk   = J[1:N,:]
    Jm   = J[N+1:2*N,:]
    indk = Array{Array{Int}}(undef,length(price)+1)
    indm = Array{Array{Int}}(undef,length(price)+1)
  
    for j = eachindex(price)
        indk[j]  = select_comp_ind(reshape(Jk[:,j],(n_par.nm,n_par.nk,n_par.ny)),0.01)   # Indexes of coefficients sorted by their absolute size
        indm[j]  = select_comp_ind(reshape(Jm[:,j],(n_par.nm,n_par.nk,n_par.ny)),0.04)

    end
    indk[17]  = select_comp_ind(reshape(Vk,(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_value)
    indm[17]  = select_comp_ind(reshape(Vm,(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_value)
    return indk, indm, J
end
function Values(price::Vector, VkSS::Array, VmSS::Array, Γ::AbstractArray, n_par, m_par)
    price   = exp.(price)
    mcw     = price[1]
    A       = price[2]
    q       = price[3]
    RL      = price[4]
    τprog   = price[5]
    τlev    = price[6]
    H       = price[7]
    Ht      = price[8]
    π       = price[9]
    r       = price[10]
    w       = price[11]
    N       = price[12]
    profits = price[13]
    unionprofits= price[14]
    av_tax_rate= price[15]
    σ       = price[16]

    Π                  = n_par.Π .+ zeros(eltype(price),1)[1]
    PP                 = ExTransition(m_par.ρ_h,n_par.bounds_y,sqrt(σ))
    Π[1:end-1,1:end-1] = PP.*(1.0-m_par.ζ)
    _, inc, eff_int = incomes(n_par, m_par, mcw, A, q, RL, τprog, τlev, H, Ht, 
                                π, r, w, N, profits, unionprofits, av_tax_rate)
    
    # expected marginal values
    EVkPrime = reshape(VkSS,(n_par.nm,n_par.nk, n_par.ny)).+ zeros(eltype(price),1)[1]
    EVmPrime = reshape(VmSS,(n_par.nm,n_par.nk, n_par.ny)).+ zeros(eltype(price),1)[1]
    @views @inbounds begin
        for mm = 1:n_par.nm
            EVkPrime[mm,:,:] .= EVkPrime[mm,:,:]*Π'
            EVmPrime[mm,:,:] .= EVmPrime[mm,:,:]*Π'
        end
    end
    # Calculate optimal policies
    c_a_star, m_a_star, k_a_star, c_n_star, m_n_star =
        EGM_policyupdate(EVmPrime ,EVkPrime ,q, π, RL.*A, 1.0, inc,n_par,m_par, false) # policy iteration
    

    # store policies
    # Update marginal values
    Vk_up, Vm_up = updateV(EVkPrime ,c_a_star, c_n_star, m_n_star, r - 1.0, q, m_par, n_par) # update expected marginal values time t
    
    Vm_up .*= eff_int
 
    Vk_up[:].= m_par.β.* Γ* Vk_up[:]
    Vm_up[:] .= m_par.β.* Γ* Vm_up[:]

    Vk_up = log.(invmutil(Vk_up,m_par))
    Vm_up = log.(invmutil(Vm_up,m_par))
    return [Vk_up[:]; Vm_up[:]]
end


function Transition(price::Vector, VkSS::Array, VmSS::Array, n_par, m_par)
    price   = exp.(price)
    mcw     = price[1]
    A       = price[2]
    q       = price[3]
    RL      = price[4]
    τprog   = price[5]
    τlev    = price[6]
    H       = price[7]
    Ht      = price[8]
    π       = price[9]
    r       = price[10]
    w       = price[11]
    N       = price[12]
    profits = price[13]
    unionprofits= price[14]
    av_tax_rate= price[15]
    σ       = price[16]

    Π                  = n_par.Π .+ zeros(eltype(price),1)[1]
    PP                 = ExTransition(m_par.ρ_h,n_par.bounds_y,sqrt(σ))
    Π[1:end-1,1:end-1] = PP.*(1.0-m_par.ζ)
    _, inc, eff_int = incomes(n_par, m_par, mcw, A, q, RL, τprog, τlev, H, Ht, 
                                π, r, w, N, profits, unionprofits, av_tax_rate)
    
    # expected marginal values
    EVkPrime = reshape(VkSS,(n_par.nm,n_par.nk, n_par.ny)).+ zeros(eltype(price),1)[1]
    EVmPrime = reshape(VmSS,(n_par.nm,n_par.nk, n_par.ny)).+ zeros(eltype(price),1)[1]
    @views @inbounds begin
        for mm = 1:n_par.nm
            EVkPrime[mm,:,:] .= EVkPrime[mm,:,:]*Π'
            EVmPrime[mm,:,:] .= EVmPrime[mm,:,:]*Π'
        end
    end
    # Calculate optimal policies
    c_a_star, m_a_star, k_a_star, c_n_star, m_n_star =
        EGM_policyupdate(EVmPrime ,EVkPrime ,q, π, RL.*A, 1.0, inc,n_par,m_par, false) # policy iteration
 
    S_a, T_a, W_a, S_n, T_n, W_n    = MakeTransition(m_a_star,  m_n_star, k_a_star, n_par.Π, n_par)
    TransitionMat_a                 = sparse(S_a, T_a, W_a, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    TransitionMat_n                 = sparse(S_n, T_n, W_n, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    Γ                               = m_par.λ .* TransitionMat_a .+ (1.0 .- m_par.λ) .* TransitionMat_n
    Vk_up, Vm_up = updateV(EVkPrime ,c_a_star, c_n_star, m_n_star, r - 1.0, q, m_par, n_par) # update expected marginal values time t
    
    Vm_up .*= eff_int

    Vk_up = log.(invmutil(Vk_up,m_par))
    Vm_up = log.(invmutil(Vm_up,m_par))
    
    return Γ, Vk_up, Vm_up
end