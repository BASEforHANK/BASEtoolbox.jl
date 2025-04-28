"""
    distrSummaries(
        distr::AbstractArray,
        q::Real,
        x_a_star::AbstractArray,
        x_n_star::AbstractArray,
        n_par,
        net_income::AbstractArray,
        gross_income::AbstractArray,
        m_par;
    )

Compute distributional summary statistics for income and wealth, including top-10% wealth,
income, and net income shares, Gini indexes for wealth and consumption, and the standard
deviation of log labor earnings.

The function calculates various metrics to summarize the distribution of wealth,
consumption, and income based on the provided data. It includes Gini coefficients, shares of
wealth and income for the top 10%, and the standard deviation of log labor earnings.

# Arguments

  - `distr::AbstractArray`: Joint distribution over liquid and illiquid assets and income
  - `q::Real`: Price of illiquid assets
  - `x_a_star::AbstractArray`: Optimal consumption policy with capital adjustment
  - `x_n_star::AbstractArray`: Optimal consumption policy without capital adjustment
  - `n_par::NumericalParameters`
  - `net_income::AbstractArray`: Vector of (on grid-)incomes (net)
  - `gross_income::AbstractArray`: Vector of (on grid-)incomes (gross)
  - `m_par::ModelParameters`

# Returns

  - `distr_b::AbstractArray`: Distribution summary over the first dimension
  - `distr_k::AbstractArray`: Distribution summary over the second dimension
  - `distr_h::AbstractArray`: Distribution summary over the third dimension
  - `TOP10Wshare::Float64`: Top 10% wealth share
  - `TOP10Ishare::Float64`: Top 10% gross income share
  - `TOP10Inetshare::Float64`: Top 10% net income share
  - `giniwealth::Float64`: Gini coefficient for wealth
  - `giniconsumption::Float64`: Gini coefficient for consumption
  - `sdlogy::Float64`: Standard deviation of log labor earnings
"""
function distrSummaries(
    distr::AbstractArray,
    q::Real,
    x_a_star::AbstractArray,
    x_n_star::AbstractArray,
    n_par,
    net_income::AbstractArray,
    gross_income::AbstractArray,
    m_par,
)

    ## ------------------------------------------------------------------------------------
    ## Step 0: Take care of complete markets case
    ## ------------------------------------------------------------------------------------

    if typeof(n_par.model) == CompleteMarkets
        distr_h = sum(distr; dims = (1, 2))[:]
        return ones(1, 1, 1),
        ones(1, 1, 1),
        distr_h,
        eps(),
        eps(),
        eps(),
        eps(),
        eps(),
        eps()
    end

    ## ------------------------------------------------------------------------------------
    ## Step 1: Compute distributional summaries for heterogeneous agents case
    ## ------------------------------------------------------------------------------------

    ## Distributional summaries
    distr_b = sum(distr; dims = (2, 3))[:]
    distr_k = sum(distr; dims = (1, 3))[:]
    distr_h = sum(distr; dims = (1, 2))[:]

    total_wealth = Array{eltype(distr)}(undef, n_par.nk .* n_par.nb)
    for k = 1:(n_par.nk)
        for b = 1:(n_par.nb)
            total_wealth[b + (k - 1) * n_par.nb] = n_par.grid_b[b] .+ q .* n_par.grid_k[k]
        end
    end

    # Wealth shares and gini
    IX = sortperm(total_wealth)
    total_wealth = total_wealth[IX]
    total_wealth_pdf = sum(distr; dims = 3)
    total_wealth_pdf = total_wealth_pdf[IX]
    total_wealth_cdf = cumsum(total_wealth_pdf)
    total_wealth_w = total_wealth .* total_wealth_pdf # weighted
    wealthshares = cumsum(total_wealth_w) ./ sum(total_wealth_w)

    TOP10Wshare = 1.0 - mylinearinterpolate(total_wealth_cdf, wealthshares, [0.9])[1]
    giniwealth = gini(total_wealth, total_wealth_pdf)

    # Consumption distribution
    c = Array{eltype(x_a_star)}(undef, (n_par.nb, n_par.nk, n_par.nh, 2))
    distr_c = similar(c)
    aux_x = net_income[5] # adjustment for labor in GHH preferences
    c[:, :, :, 1] .= x_a_star .+ aux_x # add adjustment to consumption
    c[:, :, :, 2] .= x_n_star .+ aux_x # add adjustment to consumption
    distr_c[:, :, :, 1] .= m_par.λ .* distr
    distr_c[:, :, :, 2] .= (1 - m_par.λ) .* distr

    # Gini of goods consumption
    IX = sortperm(c[:])
    c[:] .= c[IX]
    distr_c[:] .= distr_c[IX]
    giniconsumption = gini(c, distr_c)

    # Top 10 net income share
    capital_inc = net_income[2] .+ net_income[3] .- n_par.mesh_b
    Yidio = net_income[6] .+ capital_inc
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    Y_pdf = distr[IX]
    Y_cdf = cumsum(Y_pdf)
    Y_w = Yidio .* Y_pdf
    net_incomeshares = cumsum(Y_w) ./ sum(Y_w)
    TOP10Inetshare = 1.0 .- mylinearinterpolate(Y_cdf, net_incomeshares, [0.9])[1]

    # Top 10 gross income share
    Yidio = gross_income[1] .+ capital_inc
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    Y_pdf = distr[IX]
    Y_cdf = cumsum(Y_pdf)
    Y_w = Yidio .* Y_pdf
    incomeshares = cumsum(Y_w) ./ sum(Y_w)
    TOP10Ishare = 1.0 .- mylinearinterpolate(Y_cdf, incomeshares, [0.9])[1]

    # Standard deviation of log labor earnings
    Yidio = log.(gross_income[1][:, :, 1:(end - 1)])
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    distr_aux = distr[:, :, 1:(end - 1)]
    distr_aux = distr_aux ./ sum(distr_aux[:])
    Y_pdf = distr_aux[IX]

    sdlogy = sqrt(dot(Y_pdf, Yidio .^ 2) .- dot(Y_pdf, Yidio) .^ 2)

    return distr_b,
    distr_k,
    distr_h,
    TOP10Wshare,
    TOP10Ishare,
    TOP10Inetshare,
    giniwealth,
    giniconsumption,
    sdlogy
end

"""
    gini(x, pdf)

Compute the Gini coefficient of a distribution based on the values `x` and their associated
probability density function `pdf`. This implementation is intended for a discrete
probability distribution and the concrete formula follows wikipedia.

# Arguments

  - `x::AbstractArray`: Array of values representing the distribution.
  - `pdf::AbstractArray`: PDF corresponding to the values in `x`.

# Returns

  - `gini::Float64`: The Gini coefficient of the distribution.
"""
function gini(x, pdf)
    s = 0.0
    gini = 0.0
    for i in eachindex(x)
        gini -= pdf[i] * s
        s += x[i] * pdf[i]
        gini -= pdf[i] * s
    end
    gini /= s
    gini += 1.0
    return gini
end
