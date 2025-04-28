"""
    my_integrate(f::Function, a::Number, b::Number)

Perform numerical integration of a function `f` over the interval `[a, b]` using
Gauss-Chebyshev quadrature.

The function calculates the integral by evaluating `f` at the quadrature nodes and using the
corresponding weights.

# Arguments

  - `f::Function`: The function to integrate.
  - `a::Number`: The lower bound of the integration interval.
  - `b::Number`: The upper bound of the integration interval.

# Returns

  - `I::Number`: The computed integral of `f` over the interval `[a, b]`.
"""
function my_integrate(f::Function, a::Number, b::Number)
    nodes, weights = my_qnwcheb(500, a, b)
    I = weights' * f.(nodes)
end

"""
    my_qnwcheb(n::Integer, a::Number, b::Number)

Generate quadrature nodes and weights for Gauss-Chebyshev quadrature over `[a, b]`.

This function computes `n` nodes and their corresponding weights for approximating integrals
using Gauss-Chebyshev quadrature.

# Arguments

  - `n::Integer`: Number of quadrature nodes.
  - `a::Number`: The lower bound of the integration interval.
  - `b::Number`: The upper bound of the integration interval.

# Returns

  - `nodes::Vector`: Quadrature nodes, evenly distributed based on the Chebyshev polynomial.
  - `weights::Vector`: Corresponding weights for each node.
"""
function my_qnwcheb(n::Integer, a::Number, b::Number)
    nodes = (b + a) / 2 .- (b - a) / 2 .* cos.(pi / n .* (0.5:(n - 0.5)))
    weights =
        ((b - a) / n) .* (
            cos.(pi / n .* ((1:n) .- 0.5) * (2:2:(n - 1))') *
            (-2.0 ./ ((1:2:(n - 2)) .* (3:2:n))) .+ 1
        )
    return nodes, weights
end
