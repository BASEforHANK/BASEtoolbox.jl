"""
    Ksupply(
        args_hh_prob,
        n_par,
        m_par,
        Wb,
        Wk,
        distr_guess,
        net_income,
        eff_int,
    )

Calculate the aggregate savings when households face idiosyncratic income risk by first
solving the household problem using the endogenous grid method (EGM) and then finding the
stationary distribution of idiosyncratic states.

Idiosyncratic state is tuple `(b,k,y)`, where:

 1. `b`: liquid assets,
 2. `k`: illiquid assets,
 3. `y`: labor income.

This function is used in [`find_steadystate()`](@ref) to find the stationary equilibrium, as
an input to [`Kdiff()`](@ref), and in
[`BASEforHANK.PerturbationSolution.prepare_linearization()`](@ref) to prepare the
linearization of the model.

# Arguments

  - `args_hh_prob`: Vector of arguments to the household problem
  - `RK_guess`: gross real interest rate illiquid assets
  - `n_par`, `m_par`
  - `Wb`, `Wk`: guess for marginal value functions
  - `distr_guess`: guess for stationary distribution
  - `net_income`: incomes, output of functions from the IncomesETC module
  - `eff_int`: effective interest rate, output of functions from the IncomesETC module

# Returns

  - `K`,`B`: aggregate saving in illiquid (`K`) and liquid (`B`) assets
  - `Γ`,`Γ_a`,`Γ_n`: `sparse` transition matrices (average, with [`a`] or without [`n`]
    adjustment of illiquid asset)
  - `x_a_star`,`b_a_star`,`k_a_star`,`x_n_star`,`b_n_star`: optimal policies for consumption
    [`c`], liquid [`m`] and illiquid [`k`] asset, with [`a`] or without [`n`] adjustment of
    illiquid asset
  - `Wb`,`Wk`: marginal value functions
  - `distr`: ergodic steady state of `Γ`
"""
function Ksupply(
    args_hh_prob::Vector,
    n_par,
    m_par,
    Wb::AbstractArray,
    Wk::AbstractArray,
    distr_guess::AbstractArray,
    net_income::AbstractArray,
    eff_int::AbstractArray,
)
    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Step 0: Take care of complete markets case
    ## ------------------------------------------------------------------------------------

    if typeof(n_par.model) == CompleteMarkets
        @assert @isdefined(CompMarketsCapital) "Complete Markets Model requires CompMarketsCapital function."
        rSS = (1.0 .- m_par.β) ./ m_par.β  # complete markets interest rate
        K = CompMarketsCapital(rSS, m_par)
        B = m_par.ψ * K
        Wb = ones(n_par.nb, n_par.nk, n_par.nh) .* eff_int
        Wk = ones(n_par.nb, n_par.nk, n_par.nh)
        distr = reshape((n_par.Π ^ 1000)[1, :], (n_par.nb, n_par.nk, n_par.nh))
        Γ = sparse(I(2))
        Γ_a = Γ
        Γ_n = Γ
        x_a_star = ones(n_par.nb, n_par.nk, n_par.nh)
        x_n_star = ones(n_par.nb, n_par.nk, n_par.nh)
        b_a_star = ones(n_par.nb, n_par.nk, n_par.nh)
        b_n_star = ones(n_par.nb, n_par.nk, n_par.nh)
        k_a_star = ones(n_par.nb, n_par.nk, n_par.nh)
        return K,
        B,
        Γ,
        Γ_a,
        Γ_n,
        x_a_star,
        b_a_star,
        k_a_star,
        x_n_star,
        b_n_star,
        Wb,
        Wk,
        distr
    end

    ## ------------------------------------------------------------------------------------
    ## Step 1: Preallocate variables
    ## ------------------------------------------------------------------------------------

    ## Loop variables ---------------------------------------------------------------------

    dist = 9999.0
    dist1 = dist
    dist2 = dist
    count = 0

    ## Policy and value functions ---------------------------------------------------------

    nb, nk, nh = size(Wb)

    # Policy functions on exogenous grid
    b_n_star = similar(Wb)
    b_a_star = similar(Wb)
    k_a_star = similar(Wb)
    x_a_star = similar(Wb)
    x_n_star = similar(Wb)

    # Policy functions on endogenous grid, non-adjustment case
    x_tilde_n = similar(Wb)
    b_tilde_n = similar(Wb)

    # Expected marginal value functions
    EWb = similar(Wb)
    EWk = similar(Wk)

    # New marginal value functions
    Wb_new = similar(Wb)
    Wk_new = similar(Wk)

    # Inverse marginal value functions
    iWb = invmutil(Wb, m_par)
    iWk = invmutil(Wk, m_par)
    iWb_new = similar(iWb)
    iWk_new = similar(iWk)

    # Distance between old and new inverse policy functions
    D1 = similar(Wb)
    D2 = similar(Wb)

    # Difference between expected marginal value functions of assets
    E_return_diff = similar(Wb)

    # Marginal utility of consumption
    EMU = similar(Wb)

    # Marginal utility of consumption, adjustment case
    mutil_x_a = similar(Wb)

    ## Resouce grid -----------------------------------------------------------------------

    # Asset income plus liquidation value (adjustment case)
    rental_inc = net_income[2]
    liquid_asset_inc = net_income[3]
    capital_liquidation_inc = net_income[4]

    # Exogenous resource grid for the adjustment case in EGM calculated according to eq.
    # (resources adjustment)
    R_exo_a =
        reshape(rental_inc .+ liquid_asset_inc .+ capital_liquidation_inc, (nb .* nk, nh))

    ## ------------------------------------------------------------------------------------
    ## Step 2: Loop
    ## ------------------------------------------------------------------------------------

    # Iterate over marginal values until convergence
    loop_time = time()
    while dist > n_par.ϵ && count < 10000
        count += 1

        ## Take expectations given exogenous state transition -----------------------------

        #=
        Perform matrix multiplication using BLAS and reshape the result.

        Essentially this function calculates the expected marginal value of illiquid assets
        as would EWk .= reshape(reshape(Wk, (nb .* nk, nh)) * n_par.Π', (nb, nk, nh)),
        however, the BLAS `gemm!` function allows faster computation.

        The following lines use the BLAS `gemm!` function to perform fast matrix
        multiplication on the input matrix `Wk` and the parameter matrix `n_par.Π`. The
        result is stored in `EWk`. The matrices are reshaped before and after the
        multiplication to match the required dimensions.

        Arguments:
        - `Wk`: Input matrix to be multiplied, reshaped to dimensions `(nb * nk, nh)`.
        - `n_par.Π`: Parameter matrix used for multiplication.
        - `EWk`: Output matrix to store the result, reshaped to dimensions `(nb * nk, nh)`.
        - `n`: Tuple containing the dimensions for reshaping the matrices.

        Additionally, BLAS.gemm! has inputs 1.0 and 0.0 to specify potential scaling factors
        for the input matrices (not used). The arguments 'N' and 'T' specify that the
        matrices are not transposed before multiplication.

        The final result in `EWk` is reshaped back to dimensions `(nb, nk, nh)`.
        =#

        # Calculate expected marginal value of illiquid assets using BLAS following eq.
        # (ECV1)
        if typeof(n_par.model) == TwoAsset
            BLAS.gemm!(
                'N',
                'T',
                1.0,
                reshape(Wk, (nb .* nk, nh)),
                n_par.Π,
                0.0,
                reshape(EWk, (nb .* nk, nh)),
            )
            EWk .= reshape(EWk, (nb, nk, nh))
        end

        # Repeating the same process for the expected marginal value of liquid assets
        # following eq. (ECV2)
        BLAS.gemm!(
            'N',
            'T',
            1.0,
            reshape(Wb, (nb .* nk, nh)),
            n_par.Π,
            0.0,
            reshape(EWb, (nb .* nk, nh)),
        )
        EWb .= reshape(EWb, (nb, nk, nh))

        # Applying the effective interest rate to the expected marginal value of liquid
        # assets to match the expression shown in the Envelope condition in the
        # documentation (CV2)
        EWb .*= eff_int

        ## Policy update step -------------------------------------------------------------

        # Given expected marginal values, update the policy functions
        EGM_policyupdate!(
            x_a_star,
            b_a_star,
            k_a_star,
            x_n_star,
            b_n_star,
            E_return_diff,
            EMU,
            x_tilde_n,
            b_tilde_n,
            R_exo_a,
            EWb,
            EWk,
            args_hh_prob,
            net_income,
            n_par,
            m_par,
            n_par.warn_egm,
            n_par.model,
        )

        ## Marginal value update step -------------------------------------------------------

        # Given the policy functions, update the marginal value functions
        updateW!(
            Wk_new,
            Wb_new,
            mutil_x_a,
            EWk,
            x_a_star,
            x_n_star,
            b_n_star,
            args_hh_prob,
            m_par,
            n_par,
            n_par.model,
        )

        ## Calculate distance in updates and update ----------------------------------------

        # Calculate distance of inverse marginal value functions (more conservative), and
        # update marginal value functions

        if typeof(n_par.model) == TwoAsset
            invmutil!(iWk_new, Wk_new, m_par)
        end
        invmutil!(iWb_new, Wb_new, m_par)

        if typeof(n_par.model) == TwoAsset
            D1 .= iWk_new .- iWk
        end
        D2 .= iWb_new .- iWb

        if typeof(n_par.model) == TwoAsset
            dist1 = maximum(abs, D1)
        end
        dist2 = maximum(abs, D2)

        if typeof(n_par.model) == TwoAsset
            dist = max(dist1, dist2)
        else
            dist = dist2
        end

        if typeof(n_par.model) == TwoAsset
            Wk .= Wk_new
        end
        Wb .= Wb_new
        if typeof(n_par.model) == TwoAsset
            iWk .= iWk_new
        end
        iWb .= iWb_new
    end
    loop_time = time() - loop_time

    if n_par.verbose
        @printf "EGM: Iterations %d, Distance %.2e, Time %.2f\n" count dist loop_time
    end

    # Correct Wk for one asset case
    if typeof(n_par.model) == OneAsset
        Wk .= 1.0
    end

    ## ------------------------------------------------------------------------------------
    ## Step 3: Find stationary distribution
    ## ------------------------------------------------------------------------------------

    ## Define transition matrix -----------------------------------------------------------

    # Calculate inputs for sparse transition matrix
    #=

    MakeTransition from fcn_maketransition.jl returns the following variables:
    - S_a: start index for adjustment case
    - T_a: target index for adjustment case
    - W_a: weight for adjustment case
    - S_n: start index for non-adjustment case
    - T_n: target index for non-adjustment case
    - W_n: weight for non-adjustment case

    =#
    S_a, T_a, W_a, S_n, T_n, W_n =
        MakeTransition(b_a_star, b_n_star, k_a_star, n_par.Π, n_par, n_par.model)

    # Create sparse transition matrix for adjustment case of dimensions nb * nk * nh times
    # nb * nk * nh such that Γ_a[S_a[k], T_a[k]] = W_a[k].
    Γ_a = sparse(
        S_a,
        T_a,
        W_a,
        n_par.nb * n_par.nk * n_par.nh,
        n_par.nb * n_par.nk * n_par.nh,
    )

    # Create sparse transition matrix for non-adjustment case of dimensions nb * nk * nh
    # times nb * nk * nh such that Γ_n[S_n[k], T_n[k]] = W_n[k].
    Γ_n = sparse(
        S_n,
        T_n,
        W_n,
        n_par.nb * n_par.nk * n_par.nh,
        n_par.nb * n_par.nk * n_par.nh,
    )

    # Joint, probability-weighted transition matrix
    Γ = m_par.λ .* Γ_a .+ (1.0 .- m_par.λ) .* Γ_n

    ## Stationary distribution ------------------------------------------------------------

    # Calculate left-hand unit eigenvector of Γ' using eigsolve from KrylovKit
    aux = real.(eigsolve(Γ', distr_guess[:], 1)[2][1])

    # Normalize and reshape to stationary distribution (nb x nk x nh)
    distr = (reshape((aux[:]) ./ sum((aux[:])), (n_par.nb, n_par.nk, n_par.nh)))

    ## ------------------------------------------------------------------------------------
    ## Step 4: Calculate aggregate savings
    ## ------------------------------------------------------------------------------------

    # Aggregate savings in illiquid assets
    supply_illiquid = dot(sum(distr; dims = (1, 3)), n_par.grid_k)

    # Aggregate savings in liquid assets
    supply_liquid = dot(sum(distr; dims = (2, 3)), n_par.grid_b)

    if typeof(n_par.model) == TwoAsset
        K = supply_illiquid
        B = supply_liquid
    elseif typeof(n_par.model) == OneAsset
        K = (1.0 - m_par.ψ) * supply_liquid
        B = m_par.ψ * supply_liquid
    end

    ## ------------------------------------------------------------------------------------
    ## Return results
    ## ------------------------------------------------------------------------------------

    return K,
    B,
    Γ,
    Γ_a,
    Γ_n,
    x_a_star,
    b_a_star,
    k_a_star,
    x_n_star,
    b_n_star,
    Wb,
    Wk,
    distr
end
