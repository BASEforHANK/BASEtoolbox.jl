"""
    incomes(n_par, m_par, args_hh_prob)

Compute various types of net and gross income and the effective rate on liquid assets for
all agents in the economy, given all relevant arguments to the household problem (such as
prices, aggregates, taxes etc.), stored in the vector `args_hh_prob`.

Each element of the returned vector `net_inc`, `gross_inc`, and `RRi` is a 3D array with
dimensions nb x nk x nh, containing the net or gross income type or the effective rate for
an agent on the corresponding grid points (liquid assets, illiquid assets, productivity).

# Arguments

  - `n_par`, `m_par`: Parameters
  - `args_hh_prob`: Vector of arguments to the household problem

# Returns

  - `net_inc`:

     1. Net labor and union income for workers, adjusted for the composite good / net
        entrepreneurial profits for entrepreneurs
     2. Rental income from illiquid assets
     3. Liquid asset income
     4. Liquidation income from illiquid assets
     5. Transformation of composite to consumption for workers
     6. Like type 1, but without adjustment for the composite good

  - `gross_inc`:

     1. Gross labor and union income for workers / gross entrepreneurial profits for
        entrepreneurs
     2. Rental income from illiquid assets
     3. Liquid asset income
     4. Liquidation income from illiquid assets
     5. Gross labor income of workers, without union profits / gross entrepreneurial profits
        for entrepreneurs
  - `RRi`: Effective real gross return on liquid assets
"""
function incomes(n_par, m_par, args_hh_prob)

    # Number of income types
    net_inc_types = 6
    gross_inc_types = 5

    # Initialize arrays
    net_inc = fill(Array{typeof(args_hh_prob[1])}(undef, size(n_par.mesh_b)), net_inc_types)
    gross_inc =
        fill(Array{typeof(args_hh_prob[1])}(undef, size(n_par.mesh_b)), gross_inc_types)
    RRi = Array{typeof(args_hh_prob[1])}(undef, size(n_par.mesh_b))

    # Call in-place version
    incomes!(net_inc, gross_inc, RRi, n_par, m_par, args_hh_prob)

    return net_inc, gross_inc, RRi
end

"""
    incomes!(net_inc, gross_inc, RRi, n_par, m_par, args_hh_prob)

In-place version of [`incomes`](@ref), see that function for details.
"""
function incomes!(net_inc, gross_inc, RRi, n_par, m_par, args_hh_prob)
    @read_args_hh_prob()

    # Compute aggregate labor compensation
    labor_compensation = wH .* N ./ Hprog

    # Compute aggregate tax base of labor income tax
    tax_base = labor_compensation .+ Π_E

    #=

    The net income array consists of the following four crucial elements:
    1. Net labor and union income for workers, adjusted for the composite good / net
       entrepreneurial profits for entrepreneurs
    2. Rental income from illiquid assets
    3. Liquid asset income
    4. Liquidation income from illiquid assets

    These are "simply" the elements of the budget constaint in eq. (BC with x).

    Additionally, the net income array contains the following two elements:
    5. Transformation of composite to consumption for workers
    6. Like type 1, but without adjustment for the composite good

    The gross income array consists of the following elements:
    1. Gross labor and union income for workers / gross entrepreneurial profits for
       entrepreneurs
    2. Rental income from illiquid assets
    3. Liquid asset income
    4. Liquidation income from illiquid assets
    5. Gross labor income of workers, without union profits / gross entrepreneurial profits
       for entrepreneurs

    =#

    # Effective rate, see eq. (Return liquid)
    RRi .= RRL .* (n_par.mesh_b .> 0.0) .+ RRD .* (n_par.mesh_b .<= 0.0)

    # Type 2: gross/net income: rental income from illiquid assets
    rental_inc = (RK .- 1.0) .* n_par.mesh_k

    # Type 3: gross/net income: liquid asset income
    liquid_asset_inc = RRi .* n_par.mesh_b

    # Type 4: gross/net income: liquidation income from illiquid assets
    liquidation_inc = q .* n_par.mesh_k

    # Type 1: gross and net labor income of workers, see eq. (Gross income) and (Tax func)
    g_labor_inc =
        labor_compensation .* (n_par.mesh_h ./ Htilde) .^ scale_Hprog((Tprog .- 1.0), m_par)
    n_labor_inc =
        labor_tax_f.(g_labor_inc, (Tlev .- 1.0), (Tprog .- 1.0), tax_base, m_par.scale_prog)

    # Type 1: net labor income of workers, adjusted for composite good, see eq. (BC with x)
    n_labor_inc_adj = scale_GHH((Tprog .- 1.0), m_par) .* n_labor_inc

    # Type 1: gross and net union profits
    g_u_profits = n_par.frac_workers .* Π_U
    n_u_profits = union_tax_f.(g_u_profits, Tbar .- 1.0)

    # Type 1: gross and net entrepreneurial profits, see eq. (Gross income) and (Tax func)
    g_e_profits = n_par.mesh_h[:, :, end] .* Π_E
    n_e_profits =
        labor_tax_f.(g_e_profits, (Tlev .- 1.0), (Tprog .- 1.0), tax_base, m_par.scale_prog)

    # Type 5: transformation of composite to consumption for workers
    comp_labor_GHH =
        (1.0 .- scale_GHH((Tprog .- 1.0), m_par)) / (1.0 .+ (Tc .- 1.0)) .* n_labor_inc

    #=

    Next, combine all net and gross income sources.

    =#

    # Combine all net income sources, adjust for entrepreneurs
    net_inc[1] = n_labor_inc_adj .+ n_u_profits
    net_inc[1][:, :, end] .= n_e_profits
    net_inc[2] = rental_inc
    net_inc[3] = liquid_asset_inc
    net_inc[4] = liquidation_inc
    net_inc[5] = comp_labor_GHH
    net_inc[5][:, :, end] .= 0.0
    net_inc[6] = n_labor_inc .+ n_u_profits
    net_inc[6][:, :, end] .= n_e_profits

    # Combine all gross income sources, adjust for entrepreneurs
    gross_inc[1] = g_labor_inc .+ g_u_profits
    gross_inc[1][:, :, end] .= g_e_profits
    gross_inc[2] = rental_inc
    gross_inc[3] = liquid_asset_inc
    gross_inc[4] = liquidation_inc
    gross_inc[5] = g_labor_inc
    gross_inc[5][:, :, end] .= g_e_profits
end
