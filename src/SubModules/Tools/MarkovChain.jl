"""
    Tauchen(rho::Float64, N::Int; sigma::Float64 = 1.0, mue::Float64 = 0.0)

Generate a discrete approximation to an AR(1) process, following Tauchen (1987).

Uses importance sampling: each bin has probability 1/N to realize.

Authored by Christian Bayer, Uni Bonn, 03.05.2010.

# Arguments

  - `rho::Float64`: autocorrelation coefficient
  - `N::Int`: number of gridpoints
  - `sigma::Float64`: long-run variance
  - `mue::Float64`: mean of the AR(1) process

# Returns

  - `grid_vec::Vector`: state vector grid of dimension `N`
  - `P::Array{Float64, 2}`: transition matrix of dimension `N x N`
  - `bounds::Vector`: bin bounds of dimension `N + 1` #, transtype::Symbol = :importance)
"""
function Tauchen(rho::Float64, N::Int; sigma::Float64 = 1.0, mue::Float64 = 0.0)
    dis = Normal()
    pr_ij(x, bound1, bound2, rho, sigma_e) =
        pdf.(dis, x) .* (
            cdf.(dis, (bound2 - rho .* x) ./ sigma_e) -
            cdf.(dis, (bound1 - rho .* x) ./ sigma_e)
        )

    # Importance Sampling, generate equi-likely bins and the corresponding bounds, and
    # replace the outer bounds by finite numbers
    grid_probs = range(0.0; stop = 1.0, length = N + 1)
    bounds = quantile.(dis, grid_probs[1:end])
    bounds[1] = bounds[2] - 1.0e2
    bounds[end] = bounds[end - 1] + 1.0e2

    # Calculate grid-centers
    grid_vec = N * (pdf.(dis, bounds[1:(end - 1)]) - pdf.(dis, bounds[2:end]))

    # Calculate short run variance
    sigma_e = sqrt(1 - rho^2)

    # Initialize transition probability matrix
    P = fill(0.0, (N, N))

    for j = 1:N
        p(x) = pr_ij(x, bounds[j], bounds[j + 1], rho, sigma_e)

        # Exploit symmetry and evaluate integral
        for i = 1:(floor(Int, (N - 1) / 2) + 1)
            P[i, j] = my_integrate(p, bounds[i], bounds[i + 1])
        end
    end

    # Exploit symmetry
    P[(floor(Int, (N - 1) / 2) + 2):N, :] = P[(ceil(Int, (N - 1) / 2):-1:1), end:-1:1]

    # Make sure P is a probability matrix
    P = P ./ sum(P; dims = 2)

    grid_vec = grid_vec .* sigma .+ mue
    lmul!(sigma, bounds)

    return grid_vec, P, bounds
end

"""
    ExTransition(rho::Number, bounds::Array{Float64, 1}, riskscale::Number)

Calculates the transition probability matrix for a discretized state space model, similar to
the Tauchen method, using importance sampling and numerical integration with Gauss-Chebyshev
nodes.

# Arguments

  - `rho::Number`: The autoregressive parameter of the process.
  - `bounds::Array{Float64, 1}`: A 1D array specifying the grid boundaries of the state
    space.
  - `riskscale::Number`: A scaling factor that adjusts the short-run risk (variance) in the
    model.

# Returns

  - `P::Array{Float64, 2}`: The transition probability matrix, where `P[i, j]` represents
    the probability of transitioning from state `i` to state `j`.

# Notes

  - The grid `bounds` must be sorted in ascending order.    #similar to TAUCHEN
"""
function ExTransition(rho::Number, bounds::Array{Float64,1}, riskscale::Number)
    N = length(bounds) - 1

    # Calculate short run variance
    sigma_e = riskscale * sqrt(1 - rho^2)

    # Initialize transition probability matrix
    P = zeros(typeof(riskscale), N, N)

    for i = 1:(floor(Int, (N - 1) / 2) + 1)
        nodes, weights = my_qnwcheb(500, bounds[i], bounds[i + 1])

        # Exploit symmetry and evaluate integral
        for j = 1:N
            p(x) = pr_ij(x, bounds[j], bounds[j + 1], rho, sigma_e)
            P[i, j] = dot(weights, p.(nodes))
        end
    end

    # Exploit symmetry
    P[(floor(Int, (N - 1) / 2) + 2):N, :] = P[(ceil(Int, (N - 1) / 2):-1:1), end:-1:1]

    # Make sure P is a probability matrix
    P = P ./ sum(P; dims = 2)

    return P
end

"""
    pr_ij(x, bound1, bound2, rho, sigma_e)

Computes the transition probability of moving into the interval `[bound1, bound2]` starting
from the state value `x`, using a normal distribution with autoregressive parameter `rho`
and volatility `sigma_e`.

# Arguments

  - `x`: The current state value.
  - `bound1`: The lower bound of the target interval.
  - `bound2`: The upper bound of the target interval.
  - `rho`: The autoregressive parameter of the process.
  - `sigma_e`: The standard deviation (volatility) used to scale the bounds.

# Returns

  - `p`: The computed transition probability of moving into the interval `[bound1, bound2]`
    starting from the state `x`.
"""
function pr_ij(x, bound1, bound2, rho, sigma_e)
    mycdf(x) = 0.5 + 0.5 * erf.(x / sqrt(2.0))
    mypdf(x) = 1 / sqrt(2 * π) .* exp.(-x .^ 2 / 2.0)
    p =
        mypdf.(x) .*
        (mycdf.((bound2 - rho .* x) ./ sigma_e) - mycdf.((bound1 - rho .* x) ./ sigma_e))
    return p
end
