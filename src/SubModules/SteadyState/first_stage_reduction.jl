@doc raw"""
first_stage_reduction(VkSS, VmSS, TransitionMat_aSS, TransitionMat_nSS, NSS, 
    m_a_starSS, k_a_starSS, m_n_starSS, price, n_par, m_par)

Selects the DCT coefficients to be perturbed (first stage reduction, see Appendix C).

The selection is based both on the shape of the steady-state marginal values 
(in terms of log-inverse marginal utilities), and of the approximate factor representation
based on equation (C.8).

`price' is the vector of prices that appear in the household problem. 
Prices with known identical effects to other prices are left out.
`VkSS` and `VmSS` are the steady state marginal value functions, 
`TransitionMat_aSS` and `TransitionMat_nSS` are the transition matrix induced by 
the steady state policies conditional on adjustment (_a) and non-adjustment (_n).
`NSS` is the level of employment in steady state.
`m_a_starSS', `k_a_starSS', and `m_n_starSS` are the liquid and illiquid savings 
policies in steady state conditional on adjustment (_a) or non-adjustment (-n).
`n_par' and `m_par' the numerical and model parameters. 
WARNING: the order of prices is hard-coded in [`VFI()`](@ref)) 
"""
function first_stage_reduction(
    VkSS::Array,
    VmSS::Array,
    TransitionMat_aSS::SparseMatrixCSC,
    TransitionMat_nSS::SparseMatrixCSC,
    NSS::Float64,
    m_a_starSS::Array,
    k_a_starSS::Array,
    m_n_starSS::Array,
    price::Vector,
    n_par,
    m_par,
)
    #-------------------------------------------------------------------------------------
    # Calculate the Jacobian of the marginal value functions w.r.t. contemporaneous prices
    #-------------------------------------------------------------------------------------
    phi = 0.999 # approximate autocorrelation of prices

    J = ForwardDiff.jacobian(p -> VFI(p, VkSS, VmSS, NSS, n_par, m_par), log.(price))
    Numel = n_par.nm * n_par.nk * n_par.ny # Total grid size

    # Derivatives of policy functions as central finite differences
    Dk_ak = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(k_a_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Dm_ak = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(m_a_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Dk_am = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(k_a_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_m,
            1,
        )[:],
    )
    Dm_am = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(m_a_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_m,
            1,
        )[:],
    )
    Dk_nk = m_par.β .* (1.0 .- m_par.λ)
    Dm_nk = spdiagm(
        m_par.β .* (1.0 .- m_par.λ) .* centralderiv(
            reshape(m_n_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Dk_nm = 0.0
    Dm_nm = spdiagm(
        m_par.β .* (1.0 .- m_par.λ) .* centralderiv(
            reshape(m_n_starSS, (n_par.nm, n_par.nk, n_par.ny)),
            n_par.mesh_m,
            1,
        )[:],
    )

    # Joint transition matrix taking policy function marginals into account
    GammaTilde =
        [
            Dk_ak*TransitionMat_aSS Dm_ak*TransitionMat_aSS
            Dk_am*TransitionMat_aSS Dm_am*TransitionMat_aSS
        ] + [
            Dk_nk*TransitionMat_nSS Dm_nk*TransitionMat_nSS
            Dk_nm*TransitionMat_nSS Dm_nm*TransitionMat_nSS
        ]

    #-----------------------------------------------------------------------------------
    # Cummulate sum of derivatives
    #-----------------------------------------------------------------------------------
    W = (I - phi * GammaTilde) \ J
    Wk = W[1:Numel, :]         # Jacobian of the marginal value of illiquid assets
    Wm = W[(Numel+1):(2*Numel), :] # Jacobian of the marginal value of liquid assets

    # Take into account that perturbed value is written in log-inverse mutils
    TransformV(V) = log.(invmutil(V, m_par))
    Outerderivative(x) = ForwardDiff.derivative(V -> TransformV(V), x)

    CBarHat_k = Outerderivative.(VkSS[:]) .* Wk  # Chain rule
    CBarHat_m = Outerderivative.(VmSS[:]) .* Wm  # Chain rule

    #---------------------------------------------------------------------
    # Find the indexes that allow the DCTs to fit CBarHat, Vm, and Vk well 
    #---------------------------------------------------------------------
    indk = Array{Array{Int}}(undef, 2)
    indm = Array{Array{Int}}(undef, 2)

    # Calculate average absolute derivative (in DCT terms)
    Theta_m = similar(CBarHat_m)
    Theta_k = similar(CBarHat_k)
    for j in eachindex(price)
        Theta_m[:, j] = dct(reshape(CBarHat_m[:, j], (n_par.nm, n_par.nk, n_par.ny)))[:]
        Theta_k[:, j] = dct(reshape(CBarHat_k[:, j], (n_par.nm, n_par.nk, n_par.ny)))[:]
    end
    theta_m = (sum(abs.(Theta_m); dims = 2))
    theta_k = (sum(abs.(Theta_k); dims = 2))

    # Find those DCT indexes that explain the average derivative well
    indm[1] = select_ind(
        reshape(theta_m, (n_par.nm, n_par.nk, n_par.ny)),
        n_par.reduc_marginal_value,
    )
    indk[1] = select_ind(
        reshape(theta_k, (n_par.nm, n_par.nk, n_par.ny)),
        n_par.reduc_marginal_value,
    )

    # Add the indexes that fit the shape of the marginal value functions themselves well
    indk[end] = select_ind(
        dct(reshape(log.(invmutil(VkSS, m_par)), (n_par.nm, n_par.nk, n_par.ny))),
        n_par.reduc_value,
    )
    indm[end] = select_ind(
        dct(reshape(log.(invmutil(VmSS, m_par)), (n_par.nm, n_par.nk, n_par.ny))),
        n_par.reduc_value,
    )

    return indk, indm, J
end

function VFI(price::Vector, VkSS::Array, VmSS::Array, NSS::Float64, n_par, m_par)
    # Update contemporaneous value functions for given contemporaneous prices
    # and steady-state continuation marginal values for illiquid and liquid assets

    # Prices in logs (elasticities)
    price = exp.(price)

    # Read out individual "prices" (WARNING: hard coded order)
    mcw = price[1]
    q = price[2]
    RL = price[3]
    τprog = price[4]
    τlev = price[5]
    r = price[6]
    w = price[7]

    profits = price[8]
    unionprofits = price[9]
    av_tax_rate = price[10]

    σ = 1.0 + (price[11] - 1.0) * 1500 # order of magnitude adjustment

    # Fixed inputs / prices that have the same impact on decisions/Value 
    # Functions as some of the above
    A = 1.0     #  shows up only multiplicatively with RL
    π = 1.0     #  shows up only multiplicatively with RL
    H = n_par.H # in StSt: shows up only multiplicatively with w
    Ht = 1.0     # shows up only multiplicatively with w
    N = NSS     # shows up only multiplicatively with w

    # Human capital transition
    Π = n_par.Π .+ zeros(eltype(price), 1)[1]
    PP = ExTransition(m_par.ρ_h, n_par.bounds_y, sqrt(σ))
    Π[1:(end-1), 1:(end-1)] = PP .* (1.0 - m_par.ζ)

    # Incomes given prices
    _, inc, eff_int = incomes(
        n_par,
        m_par,
        mcw,
        A,
        q,
        RL,
        τprog,
        τlev,
        H,
        Ht,
        π,
        r,
        w,
        N,
        profits,
        unionprofits,
        av_tax_rate,
    )

    # Expected marginal values
    EVkPrime = reshape(VkSS, (n_par.nm, n_par.nk, n_par.ny)) .+ zeros(eltype(price), 1)[1]
    EVmPrime = reshape(VmSS, (n_par.nm, n_par.nk, n_par.ny)) .+ zeros(eltype(price), 1)[1]
    @views @inbounds begin
        for mm = 1:(n_par.nm)
            EVkPrime[mm, :, :] .= EVkPrime[mm, :, :] * Π'
            EVmPrime[mm, :, :] .= EVmPrime[mm, :, :] * Π'
        end
    end

    # Calculate optimal policies
    c_a_star, _, _, c_n_star, m_n_star =
        EGM_policyupdate(EVmPrime, EVkPrime, q, π, RL .* A, 1.0, inc, n_par, m_par, false) # policy iteration

    # Update marginal values
    Vk_up, Vm_up = updateV(EVkPrime, c_a_star, c_n_star, m_n_star, r - 1.0, q, m_par, n_par) # update expected marginal values time t

    Vm_up .*= eff_int

    return [Vk_up[:]; Vm_up[:]]
end
