"""
    EGM_policyupdate(
        EWbPrime::Array,
        EWkPrime::Array,
        args_hh_prob::Vector,
        net_income::Array,
        n_par,
        m_par,
        warnme::Bool,
        model::Union{CompleteMarkets,OneAsset,TwoAsset},
    )

Find optimal policies, given tomorrows marginal continuation values `EWbPrime`, `EWkPrime`,
today's prices [`args_hh_prob`], and income [`net_income`], using the Endogenous Grid
Method.

Optimal policies are defined over the exogenously fixed grid, while values of optimal
policies (`b` and `k`) can have off-grid values. Please refer to the subsection with the
title 'Update the optimal policy functions' of the document ['Computational
Notes.md'](Computational Notes.md), for a detailed explanation of the function's code. The
FOC's mentioned in the code are the Euler Equations as in the ['documentation of the
household problem'](HouseholdProblem.md).

# Arguments

  - `EWbPrime`, `EWkPrime`: Marginal continuation values of liquid and illiquid assets
  - `args_hh_prob`: Vector of arguments to the household problem
  - `net_income`: Incomes, output of functions from the IncomesETC module
  - `n_par`, `m_par`
  - `warnme`: If `true`, warns about non-monotonicity in resources or liquid asset choices
  - `model`: Model type, either `CompleteMarkets`, `OneAsset`, or `TwoAsset`

# Returns

  - `x_a_star`,`b_a_star`,`k_a_star`,`x_n_star`,`b_n_star`: Optimal (on-grid) policies for
    the composite [`x`], liquid [`b`] and illiquid [`k`] asset, with [`a`] or without [`n`]
    adjustment of illiquid asset
"""
function EGM_policyupdate(
    EWbPrime::Array,
    EWkPrime::Array,
    args_hh_prob::Vector,
    net_income::Array,
    n_par,
    m_par,
    warnme::Bool,
    model::Union{CompleteMarkets,OneAsset,TwoAsset},
)

    ## ------------------------------------------------------------------------------------
    ## Preallocate variables
    ## ------------------------------------------------------------------------------------

    ## Policy and value functions ---------------------------------------------------------

    nb, nk, nh = size(EWbPrime)

    # Policy functions
    b_n_star = similar(EWbPrime)
    b_a_star = similar(EWbPrime)
    k_a_star = similar(EWbPrime)
    x_a_star = similar(EWbPrime)
    x_n_star = similar(EWbPrime)

    # Policy functions on endogenous grid, non-adjustment case
    x_tilde_n = similar(EWbPrime)
    b_tilde_n = similar(EWbPrime)

    # Difference between expected marginal value functions of assets
    E_return_diff = similar(EWbPrime)

    # Marginal utility of the composite
    EMU = similar(EWbPrime)

    ## Resource grid ----------------------------------------------------------------------

    # Asset income plus liquidation value (adjustment case)
    rental_inc = net_income[2]
    liquid_asset_inc = net_income[3]
    capital_liquidation_inc = net_income[4]

    # Exogenous grid of non-human resources, calculated based on the equation (resources
    # adjustment) in the documentation
    R_exo_a =
        reshape(rental_inc .+ liquid_asset_inc .+ capital_liquidation_inc, (nb .* nk, nh))

    ## ------------------------------------------------------------------------------------
    ## Call the inplace function to update policies
    ## ------------------------------------------------------------------------------------

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
        EWbPrime,
        EWkPrime,
        args_hh_prob,
        net_income,
        n_par,
        m_par,
        warnme,
        model,
    )

    ## ------------------------------------------------------------------------------------
    ## Return results
    ## ------------------------------------------------------------------------------------

    return x_a_star, b_a_star, k_a_star, x_n_star, b_n_star
end

"""
    EGM_policyupdate!(
        x_a_star::Array,
        b_a_star::Array,
        k_a_star::Array,
        x_n_star::Array,
        b_n_star::Array,
        E_return_diff::Array,
        EMU::Array,
        x_tilde_n::Array,
        b_tilde_n::Array,
        EWbPrime::Array,
        EWkPrime::Array,
        args_hh_prob::Vector,
        net_income::Array,
        n_par,
        m_par,
        warnme::Bool,
        model::Union{CompleteMarkets,OneAsset,TwoAsset},
    )

In-place version of [`EGM_policyupdate`](@ref), see that function for details.
"""
function EGM_policyupdate!(
    x_a_star::Array,
    b_a_star::Array,
    k_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    E_return_diff::Array,
    EMU::Array,
    x_tilde_n::Array,
    b_tilde_n::Array,
    R_exo_a::Array,
    EWbPrime::Array,
    EWkPrime::Array,
    args_hh_prob::Vector,
    net_income::Array,
    n_par,
    m_par,
    warnme::Bool,
    model::CompleteMarkets,
)
    # Do nothing, complete markets does not require policy update
end

function EGM_policyupdate!(
    x_a_star::Array,
    b_a_star::Array,
    k_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    E_return_diff::Array,
    EMU::Array,
    x_tilde_n::Array,
    b_tilde_n::Array,
    R_exo_a::Array,
    EWbPrime::Array,
    EWkPrime::Array,
    args_hh_prob::Vector,
    net_income::Array,
    n_par,
    m_par,
    warnme::Bool,
    model::OneAsset,
)
    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Unpack
    ## ------------------------------------------------------------------------------------

    # Dimensions and discount factor
    nb, nk, nh = size(EWbPrime)
    β::Float64 = m_par.β

    # Income components
    net_labor_union_inc_GHH = net_income[1]
    rental_inc = net_income[2]
    liquid_asset_inc = net_income[3]

    ## ------------------------------------------------------------------------------------
    ## Part a): Find optimal liquid asset policy in the non-adjustment case
    ## ------------------------------------------------------------------------------------

    #=

    Step i): Calculate the right-hand-side of (FOC3), $EMU$, from marginal value function

    Step ii): Apply the inverse of the marginal utility function to the right-hand-side of
    (FOC3) to update the prior guess for optimal composite choice given assets tomorrow

    Step iii): Find today's endogenous grid points associated with the optimal composite
    policy choice in two steps. First, obtain assets from the budget constraint and then
    apply correct interest rate depending on borrowing or saving.

    Step iv): Interpolate back from endogenous grid onto the fixed exogenous grid

    Step v): Ensure that the borrowing constraint is respected

    =#

    # Step i): right-hand-side of Euler equation
    EMU .= EWbPrime .* β .* (1.0 .+ (Tc .- 1.0))

    # Step ii): optimal composite choice
    invmutil!(x_tilde_n, EMU, m_par)

    # Step iii): Endogenous grid points for liquid assets, calculated based on the equation
    # (end. grid non-adj.) in the documentation
    b_tilde_n .= (
        (1.0 .+ (Tc .- 1.0)) .* x_tilde_n .+ n_par.mesh_b .- net_labor_union_inc_GHH .-
        rental_inc
    )
    eff_int = RRL .* (b_tilde_n .> 0.0) + RRD .* (b_tilde_n .<= 0.0)
    b_tilde_n .= b_tilde_n ./ eff_int

    # Check monotonicity of b_tilde_n
    if warnme
        b_star_aux = reshape(b_tilde_n, (nb, nk * nh))
        if any(any(diff(b_star_aux; dims = 1) .< 0))
            @warn "Non monotone future liquid asset choice encountered."
        end
    end

    # Step iv): interpolate from endogenous grid to fixed grid
    @inbounds @views begin
        for jj = 1:nh
            for kk = 1:nk

                # Generate composite and liquid asset policies on fixed grid
                mylinearinterpolate_mult2!(
                    x_n_star[:, kk, jj],
                    b_n_star[:, kk, jj],
                    b_tilde_n[:, kk, jj],
                    x_tilde_n[:, kk, jj],
                    n_par.grid_b,
                    n_par.grid_b,
                )

                # Step v): Check for binding borrowing constraints, no extrapolation from
                # grid
                bcpol = b_tilde_n[1, kk, jj]
                for bb = 1:nb
                    if n_par.grid_b[bb] .< bcpol
                        x_n_star[bb, kk, jj] =
                            (
                                net_labor_union_inc_GHH[bb, kk, jj] .+
                                rental_inc[bb, kk, jj] .+ liquid_asset_inc[bb, kk, jj] .-
                                n_par.grid_b[1]
                            ) ./ (1 .+ (Tc .- 1.0))
                        b_n_star[bb, kk, jj] = n_par.grid_b[1]
                    end
                end
            end
        end
    end

    x_a_star .= 0.0
    b_a_star .= 0.0
    k_a_star .= 0.0
    E_return_diff .= 0.0
end

function EGM_policyupdate!(
    x_a_star::Array,
    b_a_star::Array,
    k_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    E_return_diff::Array,
    EMU::Array,
    x_tilde_n::Array,
    b_tilde_n::Array,
    R_exo_a::Array,
    EWbPrime::Array,
    EWkPrime::Array,
    args_hh_prob::Vector,
    net_income::Array,
    n_par,
    m_par,
    warnme::Bool,
    model::TwoAsset,
)
    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Unpack
    ## ------------------------------------------------------------------------------------

    # Dimensions and discount factor
    nb, nk, nh = size(EWbPrime)
    β::Float64 = m_par.β

    # Income components
    net_labor_union_inc_GHH = net_income[1]
    rental_inc = net_income[2]
    liquid_asset_inc = net_income[3]
    capital_liquidation_inc = net_income[4]

    ## ------------------------------------------------------------------------------------
    ## Part a): Find optimal liquid asset policy in the non-adjustment case
    ## ------------------------------------------------------------------------------------

    #=

    Step i): Calculate the right-hand-side of (FOC3), $EMU$, from marginal value function

    Step ii): Apply the inverse of the marginal utility function to the right-hand-side of
    (FOC3) to update the prior guess for optimal composite choice given assets tomorrow

    Step iii): Find today's endogenous grid points associated with the optimal composite
    policy choice in two steps. First, obtain assets from the budget constraint and then
    apply correct interest rate depending on borrowing or saving.

    Step iv): Interpolate back from endogenous grid onto the fixed exogenous grid

    Step v): Ensure that the borrowing constraint is respected

    =#

    # Step i): right-hand-side of Euler equation
    EMU .= EWbPrime .* β .* (1.0 .+ (Tc .- 1.0))

    # Step ii): optimal composite choice
    invmutil!(x_tilde_n, EMU, m_par)

    # Step iii): Endogenous grid points for liquid assets, calculated based on the equation
    # (end. grid non-adj.) in the documentation
    b_tilde_n .= (
        (1.0 .+ (Tc .- 1.0)) .* x_tilde_n .+ n_par.mesh_b .- net_labor_union_inc_GHH .-
        rental_inc
    )
    eff_int = RRL .* (b_tilde_n .> 0.0) + RRD .* (b_tilde_n .<= 0.0)
    b_tilde_n .= b_tilde_n ./ eff_int

    # Check monotonicity of b_tilde_n
    if warnme
        b_star_aux = reshape(b_tilde_n, (nb, nk * nh))
        if any(any(diff(b_star_aux; dims = 1) .< 0))
            @warn "Non monotone future liquid asset choice encountered."
        end
    end

    # Step iv): interpolate from endogenous grid to fixed grid
    @inbounds @views begin
        for jj = 1:nh
            for kk = 1:nk

                # Generate composite and liquid asset policies on fixed grid
                mylinearinterpolate_mult2!(
                    x_n_star[:, kk, jj],
                    b_n_star[:, kk, jj],
                    b_tilde_n[:, kk, jj],
                    x_tilde_n[:, kk, jj],
                    n_par.grid_b,
                    n_par.grid_b,
                )

                # Step v): Check for binding borrowing constraints, no extrapolation from
                # grid
                bcpol = b_tilde_n[1, kk, jj]
                for bb = 1:nb
                    if n_par.grid_b[bb] .< bcpol
                        x_n_star[bb, kk, jj] =
                            (
                                net_labor_union_inc_GHH[bb, kk, jj] .+
                                rental_inc[bb, kk, jj] .+ liquid_asset_inc[bb, kk, jj] .-
                                n_par.grid_b[1]
                            ) ./ (1 .+ (Tc .- 1.0))
                        b_n_star[bb, kk, jj] = n_par.grid_b[1]
                    end
                end
            end
        end
    end

    ## ------------------------------------------------------------------------------------
    ## Part b): Find optimal portfolio combination in the adjustment case
    ## ------------------------------------------------------------------------------------

    #=

    Step i): Calculate difference between expected marginal value functions of assets (on
    grid) on the right-hand-side of (FOC1) and (FOC2). Note that the pre-factor β * (1.0 +
    τc) does not matter for the root in the next step, hence, it is omitted.

    Step ii): Find the optimal combination of liquid asset b'*_a(k',h) and illiquid asset k'
    given a productivity state h. That is, find the b' that yields the same expected
    marginal value of liquid assets as the marginal value of illiquid assets for each
    combination of k' and h. For tractable notation we denote this as b'*(k',h) in the
    following steps of the adjustment case.

    Details on Step ii):

    - Find indifferent b by interpolation of two neighboring points m, n ∈ grid_b with
      E_return_diff(m) < 0 < E_return_diff(n)

    - Fastroot does not allow for extrapolation and uses non-negativity constraint and
      monotonicity, there is a root between the point a (defined with index idx within
      Fastroot) and the next grid point above it.

    - Fastroot returns off-grid values of b'*(k',h), but respects the borrowing constraint.

    - Example 1: Suppose E_return_diff is negative for all grid points, that means the
      expected marginal value of illiquid assets is always lower than the marginal value of
      liquid assets, therefore the household optimally chooses the highest possible liquid
      assets holding.

    - Example 2: Suppose E_return_diff is positive for all grid points, that means the
      expected marginal value of liquid assets is always lower than the marginal value of
      illiquid assets, therefore the household optimally chooses the lowest possibe liquid
      assets holding (borrowing constraint).

    =#

    # Step i) difference between expected marginal value functions of assets
    E_return_diff .= ((EWkPrime ./ q) .- EWbPrime)

    # Step ii): find optimal choice of liquid assets b' for each combination of k' and h
    # This returns an array b'*(k',h) with off-grid values
    b_hat_a1 = Fastroot(n_par.grid_b, E_return_diff)
    b_hat_a = reshape(b_hat_a1, (nk, nh))

    #=

    Step iii): Find the optimal composite policies and the endogenous grid points in the
    case that the household is only constrained in the liquid asset

    Details on Step iii):

    Ultimately in the next following steps, the optimal choices of portfolio combinations
    from above are used to derive the optimal composite policies over the endogenous grid
    points in three different cases. Step iii) deals with binding constraints for liquid
    assets and non-binding constraints for illiquid assets:

    - First (before applying inverse of marginal utility), find the expected value of the
    marginal continuation values of liquid assets at the optimal portfolio choice (off-grid)
    b' from above. For this interpolate the expression from the exogenous grid `EMU` to the
    optimal portfolio choices b'*(k',h) off-grid. The code below:

    - uses linear indexing over the dimensions (k',h)

    - assigns an index (idx) of the exogenous grid for the left value of b'*(k',h) (denoted
    by xi in the code for every iteration over (k',h)). It hereby accomodates the cases
    where b'*(k',h) is below the lowest (or above the highest) grid point.

    - calculates the distance of the optimal policy to the next grid point on the left (with
    the index idx) as left-hand-side weight s

    - uses weights s to interpolate EMU to the optimal portfolio choices b'*(k',h)

    - Then calculate composite policies that match the optimal portfolio choices using a
    standard backwards EGM step.

    - Finally, calculate the total non-human resources that are compatible with the plans of
    the household. These are the endogenous grid point for the case of adjustment.

    =#

    # Step iii): choice for when liquid asset holdings are constrained

    # Auxiliary variables to move to linear indexing
    aux_index = (0:((nk * nh) - 1)) * nb
    EMU_star = Array{eltype(b_hat_a),2}(undef, (nk, nh))
    step = diff(n_par.grid_b)

    # Interpolate EMU[b',k',h] (on-grid) over b'*(k',h) to arrive at EMU[k',h] (m-dim is
    # dropped)
    for j in eachindex(b_hat_a)

        # Find optimal choice of b' for the j-th combination of k' and h
        xi = b_hat_a[j]

        # Find indexes on grid next smallest to optimal policy and the distance s of their
        # corresponding grid values to the actual value of optimal policy xi.
        if xi .> n_par.grid_b[nb - 1]
            # If xi is above the second highest grid point, we set the left-hand-side index
            # to the second highest grid point. If xi is between the second highest and the
            # highest grid point, we use standard interpolation (see EMU_star[j] below). If
            # it is higher than the highest grid point, we use linear extrapolation with
            # negative weights on the left-hand-side, that is, s is above 1.
            idx = nb - 1
        elseif xi .<= n_par.grid_b[1]
            # If xi is below or at the lowest grid point, we set the left-hand-side index to
            # the lowest grid point. We use linear extrapolation with negative weights on
            # the right-hand-side, that is, s is equal to or below 0.
            idx = 1
        else
            # For all other interior cases, we find the next lower grid point using
            # exponential search in the locate function.
            idx = locate(xi, n_par.grid_b)
        end

        # Distance of optimal policy to next grid point as right-hand-side weight
        s = (xi .- n_par.grid_b[idx]) ./ step[idx]

        # Linear interpolation of EMU
        EMU_star[j] =
            (EMU[idx .+ aux_index[j]] .* (1.0 - s)) .+
            (s .* (EMU[idx .+ aux_index[j] .+ 1]))
    end

    # Calculate composite policies compatible with optimal portfolio choices (b'*(k',h),k')
    x_tilde_a = invmutil(EMU_star, m_par)

    # Note: Case of binding constraints for liquid assets (i.e., case where b'*(k',h) is at
    # the BC) is also considered here

    # Non-human resources that lead to assets and composite choice (endogenous grid points):
    # cash on hand of asset holdings incl. returns and liquidation value Based on eq (end.
    # grid adj.) in the documentation
    R_tilde_a =
        (1.0 .+ (Tc .- 1.0)) .* x_tilde_a .+ b_hat_a .+ capital_liquidation_inc[1, :, :] .-
        net_labor_union_inc_GHH[1, :, :]

    #=

    Step iv): Find the optimal composite policies and the endogenous grid points in the case
    that the household is constrained in the illiquid asset

    In the step above, we have obtained composite policies for the case of interior
    solutions in both assets and for the case that liquid asset holdings are constrained. We
    have also derived optimal portfolio choices for liquid asset holdings right at the
    constraint of k'=0 as b'*(k'=0,h), by conditioning on all grid points of illiquid assets
    (including zero).

    However, for some idiosyncratic states, described by h here, there can be cases where
    household are actually constrained in illiquid asset holdings. To identfy the states h,
    where there are such cases, we check whether the corresponding liquid asset choice
    b'*(k'=0,h) for state h is above the liquid asset borrowing constraint. For any lower
    resources, this implies that the household with state h would like to hold less than
    zero illiquid assets but instead can only reduces liquid asset holdings.

    The following code identifies the policies and endogenous grid point resources when
    households choose k'=0 and b'< b'*(k'=0,h) due to being constrained in illiquid assets
    as explained above, by:

    - First identifying the states h where there are such cases by checking whether
      b'*(k'=0,h) > BC

    - Making use of the equivalence of the non-adjustment case with k=0 and the constrained
      illiquid asset case under adjustment to obtain composite policies.

    - If there are no cases of constrained illiquid asset holdings for states h where the
      b'*(k'=0,h) <= BC, then the values obtained in step b) iii) using FastRoot suffice for
      these states. An empty array (to later append values from b) iii) is created for those
      states as a last step.

    =#

    # Step iv): choice for when illiquid asset holdings are constrained

    # Create lists to store composite, resources, liquid asset and illiquid assets choices
    cons_list = Array{Array{eltype(x_tilde_n)}}(undef, nh, 1)
    res_list = Array{Array{eltype(x_tilde_n)}}(undef, nh, 1)
    liq_list = Array{Array{eltype(x_tilde_n)}}(undef, nh, 1)
    cap_list = Array{Array{eltype(x_tilde_n)}}(undef, nh, 1)

    # Index for constrained liquid asset holdings, just temporary
    log_index = Vector{Bool}(undef, n_par.nb)

    # Liquid asset holdings that correspond to k' = 0 (illiquid assets constraint binds)
    b_hat_zero = b_hat_a[1, :]

    # composite of non-adjustment case that corresponds to k' = 0
    aux_x = reshape(x_tilde_n[:, 1, :], (nb, nh))

    # Income (from mesh, first and second dimension do not matter)
    aux_inc = reshape(net_labor_union_inc_GHH[1, 1, :], (1, nh))

    # Use composite at k'=0 from constrained problem, when b' is on grid, x_tilde_n is the
    # the policy function on the endogenous grid for the non-adjustment case.

    # Loop over productivity states, given that illiquid assets constraint binds
    for j = 1:nh

        # 2 cases for binding illiquid assets constraint:
        # 1. The liquid assets choice is above the borrowing constraint.
        # 2. The liquid assets choice is at or below the borrowing constraint.

        if b_hat_zero[j] > n_par.grid_b[1]

            # Case 1: We need to generate endogenous grid points here! Calculate composite
            # policies, when HHs chooses liquid asset holdings lower than b'*(k'=0,h) and
            # illiquid assets holdings k'=0 and save them in cons_list

            # Find interval of grid points below the optimal liquid assets choice. For
            # these, we need to create endogenous grid points.
            log_index .= n_par.grid_b .< b_hat_zero[j]

            # For constrained illiquid assets choices, we use the result from the
            # non-adjustment case. If I am choosing to have k' = 0 it's the same as if I was
            # not adjusting at k' = 0.
            x_k_cons = aux_x[log_index, j]
            cons_list[j] = x_k_cons # composite at k'=0 and b'< b'*(k'=0,h)

            # Required Resources: Liquid asset choice + composite - labor income resources
            # that lead to k'=0 and b' < b'*(k'=0,h)
            liquid_assets = n_par.grid_b[log_index]
            liq_list[j] = liquid_assets
            res_list[j] = liquid_assets .+ (1.0 .+ (Tc .- 1.0)) .* x_k_cons .- aux_inc[j]
            cap_list[j] = zeros(eltype(EWbPrime), sum(log_index))
        else
            # Case 2: We already covered this as part of the interior solutions in Fastroot
            # Create empty array in order to overwrite undef and prepare next step b) v)
            cons_list[j] = zeros(eltype(EWbPrime), 0)
            res_list[j] = zeros(eltype(EWbPrime), 0)
            liq_list[j] = zeros(eltype(EWbPrime), 0)
            cap_list[j] = zeros(eltype(EWbPrime), 0)
        end
    end

    #=

    Step v): Store the results (on endogenous grid points and policies) in lists

    Above, we have calculated policies in the adjustment case stored in x_tilde_a and
    b_hat_a, as well as the associated endogenous resources that are required for these
    policies.

    Moreover, we have calculated the policies and endogenous resources when we choose k'=0
    and b'<b'*(k'=0,h). In this case, we know that households already are constrained and
    have stored the policies and resources in cons_list, res_list, liq_list, and cap_list.

    Now, we want to fill the lists cons_list, res_list, liq_list, and cap_list with the
    policies for the unconstrained case, where we have calculated resources and policies in
    steps ii) and iii).

    We will use these lists to interpolate the policies from the endogenous grid defined in
    res_list back to the fixed grid defined in R_exo_a.

    =#

    # Step v): Store the results in lists

    # Reshaping policies from step ii) and step iii)
    x_tilde_a = reshape(x_tilde_a, (nk, nh))
    b_hat_a = reshape(b_hat_a, (nk, nh))

    # Fill lists with policies and resources from steps ii) and iii) by looping over
    # productivity states. This way, we ensure that the lists are filled in the correct
    # order
    for j = 1:nh
        append!(cons_list[j], x_tilde_a[:, j])
        append!(res_list[j], R_tilde_a[:, j])
        append!(liq_list[j], b_hat_a[:, j])
        append!(cap_list[j], n_par.grid_k)
    end

    ## ------------------------------------------------------------------------------------
    ## Part c): Standard backward EGM step given the optimal portfolio choice
    ## ------------------------------------------------------------------------------------

    #=

    We have obtained on-grid policies for the non-adjustment case in a). But for the the
    adjustment case, we have only obtained off-grid composite policies and portfolio
    choices. Now, we interpolate the off-grid policies back to the fixed exogenous grid for
    the case of adjustment.

    The code hereby implements the following steps:

    Step i): Interpolate back to fixed grid: Fill the policy functions (e.g. x_a_star) –
    flattened to match the dimensions of the exogenous non-human resource grid – with
    interpolated values from lists created above over the exogenous non-human resource grid
    'R_exo_a[:, j]'

    Step ii): Check for binding of both borrowing constraints (similar to last step in a)
    iv)) by:

    - First selecting all the exogenous grid points for each state h, where total resources
      are lower than the resources to only just be at the constraint in both assets (given
      by the lowest value of res_list[h]). This is done by creating a placeholder index
      `log_index2`.

    - Then for all these lower resources (where log_index2==1), the household is also
      constrained in both assets, thus portfolio choices are set to the constraint and the
      composite policy is assigned accordingly, following from the budget constraint.

    These final steps concludes the derivation of optimal policies on the exgenous liquid
    asset grid for the case of non-adjustment (x_n_star, b_n_star), and on the exogenous
    non-human resource grid for the case of adjustment (x_a_star, b_a_star, k_a_star)

    =#

    labor_inc_grid = net_labor_union_inc_GHH[1, 1, :][:]
    log_index2 = zeros(Bool, nb .* nk)

    @views @inbounds begin
        for j = 1:nh
            # Check monotonicity of resources
            if warnme
                if any(diff(res_list[j]) .< 0)
                    @warn "Non monotone resource list encountered: "
                    @printf "Sorting the ressource, consumption, liquid and capital lists. "
                    @printf "If this persists, adjust the grid! \n"
                    resort = sortperm(res_list[j])
                    res_list[j] = res_list[j][resort]
                    cons_list[j] = cons_list[j][resort]
                    liq_list[j] = liq_list[j][resort]
                    cap_list[j] = cap_list[j][resort]
                end
            end

            # Step i): Fill policy functions with interpolates of lists (from Part b, Step
            # v) with endogenous grid values) over the exogenous resource grid
            mylinearinterpolate_mult3!(
                x_a_star[:, :, j][:],
                b_a_star[:, :, j][:],
                k_a_star[:, :, j][:],
                res_list[j],
                cons_list[j],
                liq_list[j],
                cap_list[j],
                R_exo_a[:, j],
            )

            # Step ii): When at most one constraint binds: Lowest value of res_list
            # corresponds to resources needed today for optimal b'=0 and k'=0 tomorrow (&
            # according composite) Any resources on grid smaller than res_list[1] imply that
            # HHs consume all resources plus income and both constraints are binding:
            log_index2[:] .= reshape(R_exo_a[:, j], nb * nk) .< res_list[j][1]
            x_a_star[:, :, j][log_index2] .=
                (R_exo_a[log_index2, j] .+ labor_inc_grid[j] .- n_par.grid_b[1]) ./
                (1.0 .+ (Tc .- 1.0))
            b_a_star[:, :, j][log_index2] .= n_par.grid_b[1]
            k_a_star[:, :, j][log_index2] .= 0.0
        end
    end
end
