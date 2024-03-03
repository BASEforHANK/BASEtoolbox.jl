##----------------------------------------------------------------------------
## Basic Functions: Utility, marginal utility and its inverse, 
##                  Return on capital, Wages, Employment, Output, etc
##---------------------------------------------------------------------------

## Utility functions and margnal utility
function util(c::AbstractArray, m_par)
    if m_par.ξ == 1.0
        util = log.(c)
    elseif m_par.ξ == 2.0
        util = 1.0 - 1.0 ./ c
    elseif m_par.ξ == 4.0
        util = (1.0 - 1.0 ./ (c .* c .* c)) ./ 3.0
    else
        util = (c .^ (1.0 .- m_par.ξ) .- 1.0) ./ (1.0 .- m_par.ξ)
    end
    return util
end

function mutil(c::AbstractArray, m_par)
    if m_par.ξ == 1.0
        mutil = 1.0 ./ c
    elseif m_par.ξ == 2.0
        mutil = 1.0 ./ (c .^ 2)
    elseif m_par.ξ == 4.0
        mutil = 1.0 ./ ((c .^ 2) .^ 2)
    else
        mutil = c .^ (-m_par.ξ)
    end
    return mutil
end

function mutil!(mu::AbstractArray, c::AbstractArray, m_par)
    if m_par.ξ == 1.0
        mu .= 1.0 ./ c
    elseif m_par.ξ == 2.0
        mu .= 1.0 ./ (c .^ 2)
    elseif m_par.ξ == 4.0
        mu .= 1.0 ./ ((c .^ 2) .^ 2)
    else
        mu .= c .^ (-m_par.ξ)
    end
    return mu
end

function invmutil(mu, m_par)
    if m_par.ξ == 1.0
        c = 1.0 ./ mu
    elseif m_par.ξ == 2.0
        c = 1.0 ./ (sqrt.(mu))
    elseif m_par.ξ == 4.0
        c = 1.0 ./ (sqrt.(sqrt.(mu)))
    else
        c = 1.0 ./ mu .^ (1.0 ./ m_par.ξ)
    end
    return c
end

function invmutil!(c, mu, m_par)
    if m_par.ξ == 1.0
        c .= 1.0 ./ mu
    elseif m_par.ξ == 2.0
        c .= 1.0 ./ (sqrt.(mu))
    elseif m_par.ξ == 4.0
        c .= 1.0 ./ (sqrt.(sqrt.(mu)))
    else
        c .= 1.0 ./ mu .^ (1.0 ./ m_par.ξ)
    end
    return c
end

## Production functions and factor incomes
# Incomes (K:capital, Z: TFP): Interest rate = MPK.-δ, Wage = MPL, profits = Y-wL-(r+\delta)*K

output(z::Vector, m_par) = z[2] .* (z[1] .^ (m_par.α)) .* (z[3] .^ (1 - m_par.α)) #   Z .* K .^ (m_par.α) .* N .^ (1 - m_par.α)
output(K::Number, Z::Number, N::Number, m_par) = output([K; Z; N], m_par)

# Factor incomes as marginal products
interest(K::Number, Z::Number, N::Number, m_par) = Z .* m_par.α .* (K ./ N) .^ (m_par.α - 1.0) .- m_par.δ_0
wage(K::Number, Z::Number, N::Number, m_par) = Z .* (1 - m_par.α) .* (K ./ N) .^ m_par.α 

# Alternatively, using ForwardDiff instead of analytical derivatives. 
# Yields small differences in the results due to numerical precision differences. 
# interest(K::Number, Z::Number, N::Number, m_par) = ForwardDiff.gradient(z -> output(z, m_par), [K; Z; N])[1] - m_par.δ_0 # derivative instead of gradient leads to perturbation confusion
# wage(K::Number, Z::Number, N::Number, m_par) = ForwardDiff.gradient(z -> output(z, m_par), [K; Z; N])[3] # derivative instead of gradient leads to perturbation confusion

## Labor supply given wages (paid to households), resulting from leisure preferences
labor_supply(w, τlev, τprog, Ht, m_par) = ((1.0 .- τprog ) .* τlev .* w^(1.0 .- τprog)).^
                                            (1.0 ./ (m_par.γ .+ τprog)) .* Ht

labor_supply(w, m_par) = labor_supply(w, m_par.τlev, m_par.τprog, 1.0, m_par) # steady state version

## Asset markets
# price of tradable stock in steady state
qΠSS_fnc(Y::Number, RB, m_par) =
        m_par.ωΠ .* (1.0 .- 1.0 ./ m_par.μ) .* Y ./ 
        (RB ./ m_par.π .- 1 .+ m_par.ιΠ) + 1.0

# steady state payout to entrepreneurs
profitsSS_fnc(Y::Number, RB, m_par) =
        (1.0 - m_par.ωΠ) .* (1.0 .- 1.0 ./ m_par.μ) .* Y .+
        m_par.ιΠ .* (qΠSS_fnc(Y, RB, m_par) .- 1.0)

# Valuation of liquid wealth (stock)
value_liquid(B, qΠ, qΠlag) = 1.0 .+ (qΠ .- qΠlag) ./ B

##--------------------------------------------------------------------------------------
## Functions only used in the steady state calculations
##--------------------------------------------------------------------------------------
# Employment given labor supply = labor demand at a given capital stock and productivity
# This is used in the steady state calculations.
employment(K::Number, Z::Number, m_par) =
    (
        Z .* (1.0 - m_par.α) .*
        (m_par.τlev .* (1.0 - m_par.τprog )) .^ (1.0 / (1.0 - m_par.τprog )) .*
        K .^ (m_par.α)
    ) .^
    ((1.0 - m_par.τprog ) ./ (m_par.γ + m_par.τprog  + (m_par.α) .* (1 - m_par.τprog )))

# Capital intensity (K/N) given interest rate and productivity, 
# used for starting guesses in the steady state calculations
# Optimal capital intensity under Cobb Douglas
capital_intensity(r, m_par) = ((r + m_par.δ_0) ./ m_par.α .* m_par.μ)^(1.0 ./ (m_par.α .- 1))
# root finding example: capital_intensity(r, m_par) = find_zero(k -> interest(k, 1.0 / m_par.μ, 1.0, m_par) - r, 1.0)

# Capital used in production (in complete markets) at a given interest rate, taking labor supply into account
CompMarketsCapital(r, m_par) = capital_intensity(r, m_par) .* 
                               labor_supply(
                                        wage(
                                            capital_intensity(r, m_par), 1.0 ./ m_par.μ, 1.0, m_par
                                            ) ./ m_par.μw # wage markup needs to be taken into account
                                    , m_par)




        
