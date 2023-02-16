
@doc raw"""
DCT_select_V(VkSS, VmSS, TransitionMatSS, price, n_par, m_par)

Selects the DCT coefficients to be perturbed (first stage reduction, see Appendix C).

The selection is based both on the shape of the steady state marginal values 
(in terms of log-inverse marginal utilities), and of the approximate factor representation
based on equation (C.8).

`price' is the vector of prices that appear in the household problem.
`VkSS` and `VmSS` are the steady state marginal value functions, 
`TransitionMatSS` is the transition matrix induced by the steady state policies.
`n_par' and `m_par' the numerical and model parameters. 
WARNING: the order of prices is hard-coded in [`VFI()`](@ref)) 
"""
function DCT_select_V(VkSS, VmSS, TransitionMatSS, price, n_par, m_par)
    #-----------------------------------------------------------
    # Calculate the Jacobian of the marginal value functions w.r.t. contemporaneous prices
    #-----------------------------------------------------------
    
    J    = ForwardDiff.jacobian(p->VFI(p, VkSS, VmSS, n_par, m_par), log.(price))
    
    N    = n_par.nm*n_par.nk*n_par.ny # Total grid size
    Wk   = J[1:N,:]     # Jacobian of the marginal value of illiquid assets
    Wm   = J[N+1:2*N,:] # Jacobian of the marginal value of liquid assets

    #-----------------------------------------------------------
    # Cummulate sum of derivatives
    #-----------------------------------------------------------

    Jk = copy(Wk) # impact effect on Vk
    Jm = copy(Wm) # impact effect on Vm

    # Transition matrix including discounting
    trans =  m_par.β.* TransitionMatSS
    for t=1:1:10 # accumulate over T = 10 periods (c.f. C.8)
        Wk = trans * Wk
        Wm = trans * Wm
        Jk +=  Wk
        Jm +=  Wm
    end

    # Take into account that perturbed value is written in log-inverse mutils
    TransformV(V)= log.(invmutil(V,m_par))
    Outerderivative(x)  = ForwardDiff.derivative(V->TransformV(V),x)
    
    CBarHat_k = Outerderivative.(VkSS[:]) .* Jk  # Chain rule
    CBarHat_m = Outerderivative.(VmSS[:]) .* Jm  # Chain rule

    #-----------------------------------------------------------
    # Find the indexes that allow the DCTs to fit CBarHat, Vm, and Vk well 
    #-----------------------------------------------------------
    indk = Array{Array{Int}}(undef,length(price)+1)
    indm = Array{Array{Int}}(undef,length(price)+1)
  
    for j = eachindex(price)
        indk[j]  = select_comp_ind(reshape(CBarHat_k[:,j],(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_marginal_value)   
        indm[j]  = select_comp_ind(reshape(CBarHat_m[:,j],(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_marginal_value)
    end

    # Add the indexes that fit the shape of the marginal value functions themselves well
    indk[end]  = select_comp_ind(reshape(log.(invmutil(VkSS,m_par)),(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_value)
    indm[end]  = select_comp_ind(reshape(log.(invmutil(VmSS,m_par)),(n_par.nm,n_par.nk,n_par.ny)),n_par.reduc_value)

    return indk, indm, J
end

function VFI(price::Vector, VkSS::Array, VmSS::Array, n_par, m_par)
    # update contemporaneous value functions for given contemporaneous prices
    # and steady state continuation marginal values for illiquid and liquid assets
    
    # prices in logs (elasticities)
    price   = exp.(price)

    # read out individual "prices" (WARNING: hard coded order)
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

    # Human capital ztransition
    Π                  = n_par.Π .+ zeros(eltype(price),1)[1]
    PP                 = ExTransition(m_par.ρ_h,n_par.bounds_y,sqrt(σ))
    Π[1:end-1,1:end-1] = PP.*(1.0-m_par.ζ)

    # incomes given prices
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
    c_a_star, _, _, c_n_star, m_n_star =
        EGM_policyupdate(EVmPrime ,EVkPrime ,q, π, RL.*A, 1.0, inc,n_par,m_par, false) # policy iteration
    
    # Update marginal values
    Vk_up, Vm_up = updateV(EVkPrime ,c_a_star, c_n_star, m_n_star, r - 1.0, q, m_par, n_par) # update expected marginal values time t
    
    Vm_up .*= eff_int

    return [Vk_up[:]; Vm_up[:]]
end

function select_comp_ind(V,reduc_value)
    # Select the multidimensional DCT coefficients that maintain 1-reduc_value of the "energy" of V

    Theta               = dct(V)[:]                          # Discrete cosine transformation of marginal liquid asset value
    ind                 = sortperm(abs.(Theta[:]);rev=true)   # Indexes of coefficients sorted by their absolute size
    coeffs              = 1                                     # Container to store the number of retained coefficients
    # Find the important basis functions (discrete cosine) for VmSS
    while norm(Theta[ind[1:coeffs]])/norm(Theta) < 1 - reduc_value 
            coeffs      += 1                                    # add retained coefficients until only n_par.reduc_value hare of energy is lost
    end
    select_ind  = ind[1:coeffs]  
    return select_ind
end
