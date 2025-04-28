"""
    first_stage_reduction(
        WkSS,
        WbSS,
        Γ_aSS,
        Γ_nSS,
        b_a_starSS,
        k_a_starSS,
        b_n_starSS,
        args_hh_prob,
        include_list_idx,
        n_par,
        m_par,
    )

Performs the first stage reduction for the selection of DCT (Discrete Cosine Transform)
coefficients to perturb in a model. This is done as described in Appendix C of BBL.

The selection of the coefficients is based on both the shape of the steady-state marginal
value functions (in terms of log-inverse marginal utilities) and the approximate factor
representation of the model.

The function requires several steady-state inputs (such as marginal value functions, policy
functions, and transition matrices). Importantly, the aggregate inputs to the household
problem are collected in `args_hh_prob`. `include_list_idx` are the indices of the variables
that are to be perturbed; the others will remain at steady state.

# Arguments

  - `WkSS::Array`, `WbSS::Array`: Steady state marginal value functions
  - `Γ_aSS::SparseMatrixCSC`, `Γ_nSS::SparseMatrixCSC`: Transition matrices
  - `b_a_starSS::Array`, `k_a_starSS::Array`, `b_n_starSS::::Array`: Steady state policies
  - `args_hh_prob`: Vector of arguments to the household problem
  - `include_list_idx`: Vector of indices that determines which elements of ``args_hh_prob`
    are to be perturbed.
  - `m_par`, `n_par`

# Returns

  - `indk::Array{Array{Int}}`, `indb::Array{Array{Int}}`: The indices of the DCT
    coefficients for illiquid and liquid assets (`k` and `b`) that are selected for
    perturbation.
  - `J::Array`: The Jacobian of the marginal value functions with respect to contemporaneous
    inputs that are to be perturbed.
"""
function first_stage_reduction(
    WkSS::Array,
    WbSS::Array,
    Γ_aSS::SparseMatrixCSC,
    Γ_nSS::SparseMatrixCSC,
    b_a_starSS::Array,
    k_a_starSS::Array,
    b_n_starSS::Array,
    args_hh_prob,
    include_list_idx,
    n_par,
    m_par,
)
    ## ------------------------------------------------------------------------------------
    ## Step 1: Calculate the Jacobian of the marginal value functions
    ## ------------------------------------------------------------------------------------

    # Define the Jacobian of the marginal value function with respect to `inputsSS`
    Jaux = ForwardDiff.jacobian(p -> VFI(p, WkSS, WbSS, n_par, m_par), log.(args_hh_prob))

    # Scale the Jacobian of σ (hard-coded)
    σ_idx = findfirst(x -> x == "σ", args_hh_prob_names)
    Jaux[:, σ_idx] .= Jaux[:, σ_idx] .* 1500

    # Remove the co-linear variables in arguments of household problem
    J = Jaux[:, include_list_idx]

    # Total grid size
    Numel = n_par.nb * n_par.nk * n_par.nh

    ## ------------------------------------------------------------------------------------
    ## Step 2: Calculate the derivatives of policy functions as central finite differences
    ## ------------------------------------------------------------------------------------

    Dk_ak = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(k_a_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Db_ak = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(b_a_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Dk_ab = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(k_a_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_b,
            1,
        )[:],
    )
    Db_ab = spdiagm(
        m_par.β .* m_par.λ .* centralderiv(
            reshape(b_a_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_b,
            1,
        )[:],
    )
    Dk_nk = m_par.β .* (1.0 .- m_par.λ)
    Db_nk = spdiagm(
        m_par.β .* (1.0 .- m_par.λ) .* centralderiv(
            reshape(b_n_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_k,
            2,
        )[:],
    )
    Dk_nb = 0.0
    Db_nb = spdiagm(
        m_par.β .* (1.0 .- m_par.λ) .* centralderiv(
            reshape(b_n_starSS, (n_par.nb, n_par.nk, n_par.nh)),
            n_par.mesh_b,
            1,
        )[:],
    )

    ## ------------------------------------------------------------------------------------
    ## Step 3: Construct the joint transition matrix
    ## ------------------------------------------------------------------------------------

    GammaTilde =
        [
            Dk_ak*Γ_aSS Db_ak*Γ_aSS
            Dk_ab*Γ_aSS Db_ab*Γ_aSS
        ] + [
            Dk_nk*Γ_nSS Db_nk*Γ_nSS
            Dk_nb*Γ_nSS Db_nb*Γ_nSS
        ]

    ## ------------------------------------------------------------------------------------
    ## Step 4: Solve for Jacobian using the recursive relation
    ## ------------------------------------------------------------------------------------

    # Assumption: approximate autocorrelation of prices
    phi = 0.999

    # Solve the system
    W = (I - phi * GammaTilde) \ J

    # Select Jacobian of marginal value function with respect to k and b
    Wk = W[1:Numel, :]
    Wb = W[(Numel + 1):(2 * Numel), :]

    ## ------------------------------------------------------------------------------------
    ## Step 5: Transform the marginal values using the inverse marginal utility function
    ## ------------------------------------------------------------------------------------

    # Obtain outer derivative of transformation
    TransformV(V) = log.(invmutil(V, m_par))
    Outerderivative(x) = ForwardDiff.derivative(V -> TransformV(V), x)

    # Apply chain rule to compute transformed derivatives
    CBarHat_k = Outerderivative.(WkSS[:]) .* Wk
    CBarHat_m = Outerderivative.(WbSS[:]) .* Wb

    ## ------------------------------------------------------------------------------------
    ## Step 6: Perform the selection of DCT coefficients to fit the transformation well
    ## ------------------------------------------------------------------------------------

    # This step selects the coefficients that best match the average absolute derivative of
    # the marginal values in terms of their DCT representation.

    # Initialize index arrays
    indk = Array{Array{Int}}(undef, 2)
    indb = Array{Array{Int}}(undef, 2)

    # Apply the DCT to compute the transformation for CBarHat
    Theta_m = similar(CBarHat_m)
    Theta_k = similar(CBarHat_k)
    for j in axes(J, 2)
        Theta_m[:, j] = dct(reshape(CBarHat_m[:, j], (n_par.nb, n_par.nk, n_par.nh)))[:]
        Theta_k[:, j] = dct(reshape(CBarHat_k[:, j], (n_par.nb, n_par.nk, n_par.nh)))[:]
    end
    theta_m = sum(abs.(Theta_m); dims = 2)
    theta_k = sum(abs.(Theta_k); dims = 2)

    # Select DCT coefficients that explain the average derivative well
    indb[1] = select_ind(
        reshape(theta_m, (n_par.nb, n_par.nk, n_par.nh)),
        n_par.reduc_marginal_value,
    )
    indk[1] = select_ind(
        reshape(theta_k, (n_par.nb, n_par.nk, n_par.nh)),
        n_par.reduc_marginal_value,
    )

    # Add the DCT indices that match the shape of the marginal value functions well
    indk[end] = select_ind(
        dct(reshape(log.(invmutil(WkSS, m_par)), (n_par.nb, n_par.nk, n_par.nh))),
        n_par.reduc_value,
    )
    indb[end] = select_ind(
        dct(reshape(log.(invmutil(WbSS, m_par)), (n_par.nb, n_par.nk, n_par.nh))),
        n_par.reduc_value,
    )

    return indk, indb, J
end

function VFI(args_hh_prob, WkSS::Array, WbSS::Array, n_par, m_par)
    args_hh_prob = exp.(args_hh_prob)

    @read_args_hh_prob()

    # Additional definitions: borrowing rate
    RRD = borrowing_rate_ss(RRL, m_par)

    # Additional definitions: Human capital transition
    Π = n_par.Π .+ zeros(eltype(args_hh_prob), 1)[1]
    PP = ExTransition(m_par.ρ_h, n_par.bounds_h, sqrt(σ))
    Π[1:(end - 1), 1:(end - 1)] = PP .* (1.0 - m_par.ζ)

    # Net incomes of households
    @write_args_hh_prob()
    net_income, _, eff_int = incomes(n_par, m_par, args_hh_prob)

    # Expected marginal values
    EWkPrime = reshape(WkSS, (n_par.nb, n_par.nk, n_par.nh)) .+ zeros(eltype(Π), 1)[1]
    EWbPrime = reshape(WbSS, (n_par.nb, n_par.nk, n_par.nh)) .+ zeros(eltype(Π), 1)[1]

    @views @inbounds begin
        for bb = 1:(n_par.nb)
            EWkPrime[bb, :, :] .= (EWkPrime[bb, :, :] * Π')
            EWbPrime[bb, :, :] .= (EWbPrime[bb, :, :] * Π')
        end
    end

    # Calculate optimal policies
    x_a_star, _, _, x_n_star, b_n_star = EGM_policyupdate(
        EWbPrime,
        EWkPrime,
        args_hh_prob,
        net_income,
        n_par,
        m_par,
        n_par.warn_egm,
        n_par.model,
    )

    # Update marginal values
    Wk_up, Wb_up = updateW(
        EWkPrime,
        x_a_star,
        x_n_star,
        b_n_star,
        args_hh_prob,
        m_par,
        n_par,
        n_par.model,
    )

    Wb_up .*= eff_int

    return [Wk_up[:]; Wb_up[:]]
end
