@doc raw"""
    Fsys(X, XPrime, XSS, m_par, n_par, indexes, Γ, compressionIndexes, DC, IDC, DCD, IDCD)

Equilibrium error function: returns deviations from equilibrium around steady state.

Split computation into *Aggregate Part*, handled by [`Fsys_agg()`](@ref),
and *Heterogeneous Agent Part*.

# Arguments
- `X`,`XPrime`: deviations from steady state in periods t [`X`] and t+1 [`XPrime`]
- `XSS`: states and controls in steady state
- `Γ`, `DC`, `IDC`, `DCD`,`IDCD`: transformation matrices to retrieve marginal distributions [`Γ`],
    marginal value functions [`DC`,`IDC`], and the (linear) interpolant of the copula [`DCD`,`IDCD`] from deviations
- `indexes`,`compressionIndexes`: access `XSS` by variable names
    (DCT coefficients of compressed ``V_m`` and ``V_k`` in case of `compressionIndexes`)

# Example
```jldoctest
julia> # Solve for steady state, construct Γ,DC,IDC as in LinearSolution()
julia> Fsys(zeros(ntotal),zeros(ntotal),XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC)
*ntotal*-element Array{Float64,1}:
 0.0
 0.0
 ...
 0.0
```
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
    # The function call with Duals takes
    # Reserve space for error terms
    F = zeros(eltype(X), size(X))

    ############################################################################
    #            I. Read out argument values                                   #
    ############################################################################

    ############################################################################
    # I.1. Generate code that reads aggregate states/controls
    #      from steady state deviations. Equations take the form of:
    # r       = exp.(XSS[indexes.rSS] .+ X[indexes.r])
    # rPrime  = exp.(XSS[indexes.rSS] .+ XPrime[indexes.r])
    ############################################################################

    # @generate_equations(aggr_names)
    @generate_equations()

    ############################################################################
    # I.2. Read out  perturbed distributions
    ############################################################################

    # Copula parameters (deviations from steads state)
    θD = uncompress(compressionIndexes[3], X[indexes.COP], DCD, IDCD)
    COP_Dev = reshape(copy(θD[:]), (n_par.nm_copula, n_par.nk_copula, n_par.ny_copula))
    COP_Dev = pdf_to_cdf(COP_Dev)

    θDPrime = uncompress(compressionIndexes[3], XPrime[indexes.COP], DCD, IDCD)
    COP_DevPrime =
        reshape(copy(θDPrime), (n_par.nm_copula, n_par.nk_copula, n_par.ny_copula))
    COP_DevPrime = pdf_to_cdf(COP_DevPrime)

    # marginal distributions (pdfs, including steady state)
    distr_m = XSS[indexes.distr_mSS] .+ Γ[1] * X[indexes.distr_m]
    distr_k = XSS[indexes.distr_kSS] .+ Γ[2] * X[indexes.distr_k]
    distr_y = XSS[indexes.distr_ySS] .+ Γ[3] * X[indexes.distr_y]

    distr_m_Prime = XSS[indexes.distr_mSS] .+ Γ[1] * XPrime[indexes.distr_m]
    distr_k_Prime = XSS[indexes.distr_kSS] .+ Γ[2] * XPrime[indexes.distr_k]
    distr_y_Prime = XSS[indexes.distr_ySS] .+ Γ[3] * XPrime[indexes.distr_y]

    # marginal distributions (cdfs) 
    CDF_m = cumsum(distr_m[:])
    CDF_k = cumsum(distr_k[:])
    CDF_y = cumsum(distr_y[:])

    ############################################################################
    # I.3. Read out steady state distributions
    ############################################################################

    # steads state cdfs (on value grid)
    CDF_mSS = cumsum(XSS[indexes.distr_mSS]) .+ zeros(eltype(θD), n_par.nm)
    CDF_kSS = cumsum(XSS[indexes.distr_kSS]) .+ zeros(eltype(θD), n_par.nk)
    CDF_ySS = cumsum(XSS[indexes.distr_ySS]) .+ zeros(eltype(θD), n_par.ny)

    # steady state copula (on copula grid)
    COPSS =
        reshape(XSS[indexes.COPSS] .+ zeros(eltype(θD), 1), (n_par.nm, n_par.nk, n_par.ny))
    COPSS = pdf_to_cdf(COPSS)

    # steady state copula marginals (cdfs) 
    s_m_m = n_par.copula_marginal_m .+ zeros(eltype(θD), 1)
    s_m_k = n_par.copula_marginal_k .+ zeros(eltype(θD), 1)
    s_m_y = n_par.copula_marginal_y .+ zeros(eltype(θD), 1)

    ############################################################################
    # I.4. Produce perturbed joint distribution using the copula
    ############################################################################
    # Copula(x::AbstractVector,y::AbstractVector,z::AbstractVector) = 
    # myAkimaInterp3(CDF_mSS, CDF_kSS, CDF_ySS, COPSS, x, y, z) .+
    # myAkimaInterp3(s_m_m, s_m_k, s_m_y, COP_Dev, x, y, z)

    Copula(x::Vector, y::Vector, z::Vector) =
        myinterpolate3(CDF_mSS, CDF_kSS, CDF_ySS, COPSS, x, y, z) .+
        myinterpolate3(s_m_m, s_m_k, s_m_y, COP_Dev, x, y, z)

    CDF_joint = Copula(CDF_m[:], CDF_k[:], CDF_y[:]) # roughly 5% of time
    PDF_joint = cdf_to_pdf(CDF_joint)

    ############################################################################
    # I.5 uncompressing policies/value functions
    ###########################################################################
    VmSS = XSS[indexes.VmSS]
    VkSS = XSS[indexes.VkSS]
    VmPrime = mutil(
        exp.(VmSS .+ uncompress(compressionIndexes[1], XPrime[indexes.Vm], DC, IDC)),
        m_par,
    )
    VkPrime = mutil(
        exp.(VkSS .+ uncompress(compressionIndexes[2], XPrime[indexes.Vk], DC, IDC)),
        m_par,
    )

    ############################################################################
    #           II. Auxiliary Variables                                        #
    ############################################################################
    # Transition Matrix Productivity

    Π = n_par.Π .+ zeros(eltype(X), 1)[1]
    PP = ExTransition(m_par.ρ_h, n_par.bounds_y, sqrt(σ))
    Π[1:end-1, 1:end-1] = PP .* (1.0 - m_par.ζ)

    ############################################################################
    #           III. Error term calculations (i.e. model starts here)          #
    ############################################################################

    ############################################################################
    #           III. 1. Aggregate Part #
    ############################################################################
    F = Fsys_agg(X, XPrime, XSS, PDF_joint, m_par, n_par, indexes)

    # Average Human Capital =
    # average productivity (at the productivit grid, used to normalize to 0)
    tax_prog_scale  = (m_par.γ + m_par.τprog ) / ((m_par.γ + τprog))
    H               = dot(distr_y[1:end-1], n_par.grid_y[1:end-1])
    KP              = dot(n_par.grid_k, distr_k[:])
    Htact           = dot(distr_y[1:end-1], (n_par.grid_y[1:end-1] / H) .^ (tax_prog_scale))
    BP              = dot(n_par.grid_m, distr_m[:])
    BDact           = -sum(distr_m .* (n_par.grid_m .< 0) .* n_par.grid_m)

    ############################################################################
    #               III. 2. Heterogeneous Agent Part                           #
    ############################################################################
    # Incomes
    incgross, inc, eff_int = incomes(
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

    # Calculate Taxes
    tax_prog_scale  = (m_par.γ + m_par.τprog ) / ((m_par.γ + τprog))
    LC              = mcw * w .* N ./ Ht
    taxrev =
        ((n_par.grid_y / H) .^ tax_prog_scale .* LC) -
        τlev .* ((n_par.grid_y / H) .^ tax_prog_scale .* LC) .^ (1.0 - τprog)
    taxrev[end] =
        n_par.grid_y[end] .* profits -
        τlev .* (n_par.grid_y[end] .* profits) .^ (1.0 - τprog)
    incgrossaux      = ((n_par.grid_y / H) .^ tax_prog_scale .* LC)
    incgrossaux[end] = n_par.grid_y[end] .* profits
    av_tax_rate_up   = dot(distr_y, taxrev) ./ (dot(distr_y, incgrossaux))

    # Calculate optimal policies
    # expected margginal values
    EVkPrime = reshape(VkPrime, (n_par.nm, n_par.nk, n_par.ny))
    EVmPrime = reshape(VmPrime, (n_par.nm, n_par.nk, n_par.ny))

    @views @inbounds begin
        for mm = 1:n_par.nm
            EVkPrime[mm, :, :] .= EVkPrime[mm, :, :] * Π'
            EVmPrime[mm, :, :] .= (EVmPrime[mm, :, :] * Π')
        end
    end
    c_a_star, m_a_star, k_a_star, c_n_star, m_n_star =
        EGM_policyupdate(EVmPrime, EVkPrime, q, π, RL .* A, 1.0, inc, n_par, m_par, false) # policy iteration

    # Update marginal values
    Vk_err, Vm_err =
        updateV(EVkPrime, c_a_star, c_n_star, m_n_star, r - 1.0, q, m_par, n_par) # update expected marginal values time t
    Vm_err .*= eff_int
    if !only_F
        Vm_new = copy(Vm_err)
        Vk_new = copy(Vk_err)
    end
    # Update distribution
    dist_aux = DirectTransition(m_a_star, m_n_star, k_a_star, PDF_joint, m_par.λ, Π, n_par)
    PDF_jointPrime = reshape(dist_aux, n_par.nm, n_par.nk, n_par.ny)

    #----------------------------------------------------------------------------------------
    # Calculate Error Terms
    #----------------------------------------------------------------------------------------
    # Error terms on marginal values (controls)
    invmutil!(Vm_err, Vm_err, m_par)
    invmutil!(Vk_err, Vk_err, m_par)
    Vm_err .= log.(Vm_err) .- reshape(VmSS, (n_par.nm, n_par.nk, n_par.ny))
    Vm_thet = compress(compressionIndexes[1], Vm_err, DC, IDC)
    F[indexes.Vm] = X[indexes.Vm] .- Vm_thet

    Vk_err .= log.(Vk_err) .- reshape(VkSS, (n_par.nm, n_par.nk, n_par.ny))
    Vk_thet = compress(compressionIndexes[2], Vk_err, DC, IDC)
    F[indexes.Vk] = X[indexes.Vk] .- Vk_thet

    # Error Terms on marginal distribution (in levels, states)
    distr_mPrimeUpdate = dropdims(sum(PDF_jointPrime, dims = (2, 3)), dims = (2, 3))
    distr_kPrimeUpdate = dropdims(sum(PDF_jointPrime, dims = (1, 3)), dims = (1, 3))
    distr_yPrimeUpdate = (distr_y'*Π)[:]
    F[indexes.distr_m] = (distr_mPrimeUpdate.-distr_m_Prime)[1:end-1]
    F[indexes.distr_k] = (distr_kPrimeUpdate.-distr_k_Prime)[1:end-1]
    F[indexes.distr_y] = (distr_yPrimeUpdate.-distr_y_Prime[:])[1:end-1]

    # Error Terms on Copula (states)
    # Deviation of iterated copula from fixed copula
    # CopulaDevPrime(x::AbstractVector,y::AbstractVector,z::AbstractVector) = 
    # myAkimaInterp3(CDF_m_PrimeUp, CDF_k_PrimeUp, CDF_y_PrimeUp, pdf_to_cdf(PDF_jointPrime), x, y, z) .-
    # myAkimaInterp3(CDF_mSS, CDF_kSS, CDF_ySS, COPSS, x, y, z)
    CDF_m_PrimeUp = cumsum(distr_mPrimeUpdate)
    CDF_k_PrimeUp = cumsum(distr_kPrimeUpdate)
    CDF_y_PrimeUp = cumsum(distr_yPrimeUpdate)
    CopulaDevPrime(x::Vector, y::Vector, z::Vector) =
        myinterpolate3(
            CDF_m_PrimeUp,
            CDF_k_PrimeUp,
            CDF_y_PrimeUp,
            pdf_to_cdf(PDF_jointPrime),
            x,
            y,
            z,
        ) .- myinterpolate3(CDF_mSS, CDF_kSS, CDF_ySS, COPSS, x, y, z)

    CDF_Dev = CopulaDevPrime(s_m_m, s_m_k, s_m_y) # interpolate deviations on copula grid
    COP_thet =
        compress(compressionIndexes[3], cdf_to_pdf(CDF_Dev - COP_DevPrime), DCD, IDCD) # calculate DCT of deviations

    F[indexes.COP] = COP_thet

    # Calculate distribution statistics (generalized moments)
    _, _, _, TOP10WshareT, TOP10IshareT, TOP10InetshareT, GiniWT, GiniCT, sdlogyT =
        distrSummaries(PDF_joint, q, c_a_star, c_n_star, n_par, inc, incgross, m_par)

    # Error Term on prices/aggregate summary vars (logarithmic, controls)
    F[indexes.K]  = log.(K) - log.(KP)
    F[indexes.B]  = log.(B) - log.(BP)
    F[indexes.BD] = log.(BD) - log.(BDact)
    F[indexes.Ht] = log.(Ht) - log.(Htact)
    F[indexes.τlev] = av_tax_rate - av_tax_rate_up
    F[indexes.T]  = log(T) - log(dot(distr_y, taxrev) + av_tax_rate * (unionprofits))

    # Error Terms on  distribution summaries
    F[indexes.GiniW] = log.(GiniW) - log.(GiniWT)
    F[indexes.TOP10Ishare] = log.(TOP10Ishare) - log.(TOP10IshareT)
    F[indexes.TOP10Inetshare] = log.(TOP10Inetshare) - log.(TOP10InetshareT)
    F[indexes.TOP10Wshare] = log.(TOP10Wshare) - log.(TOP10WshareT)
    F[indexes.GiniC] = log.(GiniC) - log.(GiniCT)
    F[indexes.sdlogy] = log.(sdlogy) - log.(sdlogyT)


    if only_F
        return F
    else
        return F, c_a_star, m_a_star, k_a_star, c_n_star, m_n_star, Vk_new, Vm_new, taxrev
    end
end
