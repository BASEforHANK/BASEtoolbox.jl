"""
    Fsys(X, XPrime, XSS, m_par, n_par, indexes, Γ, compressionIndexes, DC, IDC, DCD, IDCD; only_F = true)

Equilibrium error function: returns deviations from equilibrium around steady state.

Split computation into *Aggregate Part*, handled by [`Fsys_agg()`](@ref) and *Heterogeneous
Agent Part*.

# Arguments

  - `X`,`XPrime`: deviations from steady state in periods t [`X`] and t+1 [`XPrime`]

  - `XSS`: states and controls in steady state
  - `Γ`, `DC`, `IDC`, `DCD`,`IDCD`: transformation matrices to retrieve:

      + marginal distributions [`Γ`],
      + marginal value functions [`DC`,`IDC`], and
      + the (linear) interpolant of the copula [`DCD`,`IDCD`] from deviations
  - `indexes`,`compressionIndexes`: access `XSS` by variable names (DCT coefficients of
    compressed ``V_m`` and ``V_k`` in case of `compressionIndexes`)
"""
function Fsys(
    X::AbstractArray,
    XPrime::AbstractArray,
    XSS::Array{Float64,1},
    m_par,
    n_par,
    indexes,
    Γ::Array{Array{Float64,2},1},
    compressionIndexes::Array{Array{Int,1},1},
    DC::Array{Array{Float64,2},1},
    IDC::Array{Adjoint{Float64,Array{Float64,2}},1},
    DCD::Array{Array{Float64,2},1},
    IDCD::Array{Adjoint{Float64,Array{Float64,2}},1};
    only_F = true,
)

    ## ------------------------------------------------------------------------------------
    ## Preamble
    ## ------------------------------------------------------------------------------------

    # Initialize the output vector, use the same type as the input
    F = zeros(eltype(X), size(X))

    ## Unpack aggregate variables -------------------------------------------------------

    # Unpack X, XPrime, and XSS into variables
    @generate_equations()

    ## Unpack perturbed distributions -----------------------------------------------------

    # Copula parameters (deviations from steady state)
    θD = uncompress(compressionIndexes[3], X[indexes.COP], DCD, IDCD)
    COP_Dev = reshape(copy(θD[:]), (n_par.nb_copula, n_par.nk_copula, n_par.nh_copula))
    COP_Dev = pdf_to_cdf(COP_Dev)

    θDPrime = uncompress(compressionIndexes[3], XPrime[indexes.COP], DCD, IDCD)
    COP_DevPrime =
        reshape(copy(θDPrime), (n_par.nb_copula, n_par.nk_copula, n_par.nh_copula))
    COP_DevPrime = pdf_to_cdf(COP_DevPrime)

    # marginal distributions (pdfs, including steady state)
    distr_b = XSS[indexes.distr_bSS] .+ Γ[1] * X[indexes.distr_b]
    distr_k = XSS[indexes.distr_kSS] .+ Γ[2] * X[indexes.distr_k]
    distr_h = XSS[indexes.distr_hSS] .+ Γ[3] * X[indexes.distr_h]

    distr_b_Prime = XSS[indexes.distr_bSS] .+ Γ[1] * XPrime[indexes.distr_b]
    distr_k_Prime = XSS[indexes.distr_kSS] .+ Γ[2] * XPrime[indexes.distr_k]
    distr_h_Prime = XSS[indexes.distr_hSS] .+ Γ[3] * XPrime[indexes.distr_h]

    # marginal distributions (cdfs)
    CDF_b = cumsum(distr_b[:])
    CDF_k = cumsum(distr_k[:])
    CDF_h = cumsum(distr_h[:])

    ## Unpack steady state distributions --------------------------------------------------

    # steads state cdfs (on value grid)
    CDF_bSS = cumsum(XSS[indexes.distr_bSS]) .+ zeros(eltype(θD), n_par.nb)
    CDF_kSS = cumsum(XSS[indexes.distr_kSS]) .+ zeros(eltype(θD), n_par.nk)
    CDF_hSS = cumsum(XSS[indexes.distr_hSS]) .+ zeros(eltype(θD), n_par.nh)

    # steady state copula (on copula grid)
    COPSS =
        reshape(XSS[indexes.COPSS] .+ zeros(eltype(θD), 1), (n_par.nb, n_par.nk, n_par.nh))
    COPSS = pdf_to_cdf(COPSS)

    # steady state copula marginals (cdfs)
    s_m_b = n_par.copula_marginal_b .+ zeros(eltype(θD), 1)
    s_m_k = n_par.copula_marginal_k .+ zeros(eltype(θD), 1)
    s_m_h = n_par.copula_marginal_h .+ zeros(eltype(θD), 1)

    ## Joint distribution -----------------------------------------------------------------
    Copula(x::Vector, y::Vector, z::Vector) =
        myinterpolate3(CDF_bSS, CDF_kSS, CDF_hSS, COPSS, n_par.model, x, y, z) .+
        myinterpolate3(s_m_b, s_m_k, s_m_h, COP_Dev, n_par.model, x, y, z)

    CDF_joint = Copula(CDF_b[:], CDF_k[:], CDF_h[:])
    PDF_joint = cdf_to_pdf(CDF_joint)

    ## Unpack value functions -------------------------------------------------------------

    WbSS = XSS[indexes.WbSS]
    WkSS = XSS[indexes.WkSS]
    WbPrime = mutil(
        exp.(WbSS .+ uncompress(compressionIndexes[1], XPrime[indexes.Wb], DC, IDC)),
        m_par,
    )
    WkPrime = mutil(
        exp.(WkSS .+ uncompress(compressionIndexes[2], XPrime[indexes.Wk], DC, IDC)),
        m_par,
    )

    ## Unpack transition matrix -----------------------------------------------------------

    Π = n_par.Π .+ zeros(eltype(X), 1)[1]
    if typeof(n_par.model) == OneAsset || typeof(n_par.model) == TwoAsset
        PP = ExTransition(m_par.ρ_h, n_par.bounds_h, sqrt(σ))
        Π[1:(end - 1), 1:(end - 1)] = PP .* (1.0 - m_par.ζ)
    end

    ## ------------------------------------------------------------------------------------
    ## Equilibrium conditions (aggregate)
    ## ------------------------------------------------------------------------------------

    ## Aggregate equations ----------------------------------------------------------------
    # Load aggregate equations as specified in the model file (and precompiled)
    F = Fsys_agg(X, XPrime, XSS, PDF_joint, m_par, n_par, indexes)

    ## Update distributional statistics ---------------------------------------------------
    # These are also in Fsys_agg, but we need to update them here, see documentation

    # Scaling factor for individual productivity
    F[indexes.Htilde] =
        (log(Htilde)) - (log(dot(distr_h[1:(end - 1)], n_par.grid_h[1:(end - 1)])))

    if typeof(n_par.model) == OneAsset

        # Total assets
        F[indexes.TotalAssets] = (log(TotalAssets)) - log(dot(n_par.grid_b, distr_b[:]))

        # IOUs
        F[indexes.BD] =
            (log(BD)) - (log(-sum(distr_b .* (n_par.grid_b .< 0) .* n_par.grid_b)))

    elseif typeof(n_par.model) == TwoAsset

        # Capital market clearing
        F[indexes.K] = (log(K)) - (log(dot(n_par.grid_k, distr_k[:])))

        # Bond market clearing
        F[indexes.B] = (log(B)) - (log(dot(n_par.grid_b, distr_b[:])))

        # IOUs
        F[indexes.BD] =
            (log(BD)) - (log(-sum(distr_b .* (n_par.grid_b .< 0) .* n_par.grid_b)))

    elseif typeof(n_par.model) == CompleteMarkets

        # Do nothing, everything defined in the aggregate model

    end

    ## ------------------------------------------------------------------------------------
    ## Equilibrium conditions (idiosyncratic)
    ## ------------------------------------------------------------------------------------

    ## Incomes ----------------------------------------------------------------------------
    # Calculate incomes based on the model-specific income functions

    @write_args_hh_prob()

    # Calculate net income and effective interest rate
    net_income, gross_income, eff_int = incomes(n_par, m_par, args_hh_prob)

    ## Policy and value functions ---------------------------------------------------------

    # Calculate expected marginal value functions
    EWkPrime = reshape(WkPrime, (n_par.nb, n_par.nk, n_par.nh))
    EWbPrime = reshape(WbPrime, (n_par.nb, n_par.nk, n_par.nh))
    @views @inbounds begin
        for bb = 1:(n_par.nb)
            EWkPrime[bb, :, :] .= (EWkPrime[bb, :, :] * Π')
            EWbPrime[bb, :, :] .= (EWbPrime[bb, :, :] * Π')
        end
    end

    # Calculate policy functions (policy iteration)
    x_a_star, b_a_star, k_a_star, x_n_star, b_n_star = EGM_policyupdate(
        EWbPrime,
        EWkPrime,
        args_hh_prob,
        net_income,
        n_par,
        m_par,
        false,
        n_par.model,
    )

    # Update marginal values
    Wk_err, Wb_err = updateW(
        EWkPrime,
        x_a_star,
        x_n_star,
        b_n_star,
        args_hh_prob,
        m_par,
        n_par,
        n_par.model,
    )

    Wb_err .*= eff_int
    if !only_F
        Wb_new = copy(Wb_err)
        Wk_new = copy(Wk_err)
    end

    # Update distribution
    dist_aux = DirectTransition(b_a_star, b_n_star, k_a_star, PDF_joint, m_par.λ, Π, n_par)
    PDF_jointPrime = reshape(dist_aux, n_par.nb, n_par.nk, n_par.nh)

    ## Set up the error terms for idiosyncratic part --------------------------------------

    # Error terms on marginal values (controls)
    invmutil!(Wb_err, Wb_err, m_par)
    invmutil!(Wk_err, Wk_err, m_par)
    Wb_err .= log.(Wb_err) .- reshape(WbSS, (n_par.nb, n_par.nk, n_par.nh))
    Wk_err .= log.(Wk_err) .- reshape(WkSS, (n_par.nb, n_par.nk, n_par.nh))
    Wb_thet = compress(compressionIndexes[1], Wb_err, DC, IDC)
    Wk_thet = compress(compressionIndexes[2], Wk_err, DC, IDC)
    F[indexes.Wb] = X[indexes.Wb] .- Wb_thet
    F[indexes.Wk] = X[indexes.Wk] .- Wk_thet

    # Error Terms on marginal distribution (in levels, states)
    distr_bPrimeUpdate = dropdims(sum(PDF_jointPrime; dims = (2, 3)); dims = (2, 3))
    distr_kPrimeUpdate = dropdims(sum(PDF_jointPrime; dims = (1, 3)); dims = (1, 3))
    distr_hPrimeUpdate = (distr_h' * Π)[:]
    F[indexes.distr_b] = (distr_bPrimeUpdate .- distr_b_Prime)[1:(end - 1)]
    F[indexes.distr_k] = (distr_kPrimeUpdate .- distr_k_Prime)[1:(end - 1)]
    F[indexes.distr_h] = (distr_hPrimeUpdate .- distr_h_Prime[:])[1:(end - 1)]

    # Error Terms on copula (states): deviation of iterated copula from fixed copula
    CDF_b_PrimeUp = cumsum(distr_bPrimeUpdate)
    CDF_k_PrimeUp = cumsum(distr_kPrimeUpdate)
    CDF_h_PrimeUp = cumsum(distr_hPrimeUpdate)
    CopulaDevPrime(x::Vector, y::Vector, z::Vector) =
        myinterpolate3(
            CDF_b_PrimeUp,
            CDF_k_PrimeUp,
            CDF_h_PrimeUp,
            pdf_to_cdf(PDF_jointPrime),
            n_par.model,
            x,
            y,
            z,
        ) .- myinterpolate3(CDF_bSS, CDF_kSS, CDF_hSS, COPSS, n_par.model, x, y, z)
    CDF_Dev = CopulaDevPrime(s_m_b, s_m_k, s_m_h)
    COP_thet =
        compress(compressionIndexes[3], cdf_to_pdf(CDF_Dev - COP_DevPrime), DCD, IDCD)
    F[indexes.COP] = COP_thet

    ## ------------------------------------------------------------------------------------
    ## Equilibrium conditions (auxiliary statistics)
    ## ------------------------------------------------------------------------------------

    # Calculate distribution statistics (generalized moments)
    _, _, _, TOP10WshareT, TOP10IshareT, TOP10InetshareT, GiniWT, GiniCT, sdlogyT =
        distrSummaries(
            PDF_joint,
            q,
            x_a_star,
            x_n_star,
            n_par,
            net_income,
            gross_income,
            m_par,
        )

    # Error Terms on  distribution summaries
    F[indexes.GiniW] = log.(GiniW) - log.(GiniWT)
    F[indexes.TOP10Ishare] = log.(TOP10Ishare) - log.(TOP10IshareT)
    F[indexes.TOP10Inetshare] = log.(TOP10Inetshare) - log.(TOP10InetshareT)
    F[indexes.TOP10Wshare] = log.(TOP10Wshare) - log.(TOP10WshareT)
    F[indexes.GiniC] = log.(GiniC) - log.(GiniCT)
    F[indexes.sdlogy] = log.(sdlogy) - log.(sdlogyT)

    ## ------------------------------------------------------------------------------------
    ## Return
    ## ------------------------------------------------------------------------------------

    if only_F
        return F
    else
        return F,
        x_a_star,
        b_a_star,
        k_a_star,
        x_n_star,
        b_n_star,
        Wk_new,
        Wb_new,
        ((Tbar .- 1.0) * (wH * N) + (Tbar .- 1.0) * Π_E)
    end
end
