@doc raw"""
    Brent(g::Function, a::Real, b::Real; tol = 1e-14)

Find the root of a function using Brent's method, a combination of the bisection, secant,
and inverse quadratic interpolation methods.

This implementation solves for a root of the function `g(z)`, where the function values at
`a` and `b` must have opposite signs, ensuring the existence of a root between them. The
implementation follows wikipedia.

# Arguments
- `g::Function`: The function for which the root is being found. It should return a value
    that can be indexed, where the first element (`g(z)[1]`) represents the function value
    at `z`.
- `a::Real`: The lower bound of the interval where the root is located.
- `b::Real`: The upper bound of the interval where the root is located.
- `tol::Real`: The tolerance for convergence. The default is `1e-14`.

# Returns
- `b::Real`: The estimated root of the function.
- `iter::Int`: The number of iterations used to find the root.

# Errors
- Throws an error if the function values at `a` and `b` have the same sign, indicating no
  root exists in the specified interval.

# Example
```julia
    g(x) = x^2 - 2
    root, iterations = Brent(g, 0, 2)
```
"""
function Brent(g::Function, a::Real, b::Real; tol = 1e-14)
    f(z) = g(z)[1]
    fa = f(a)
    fb = f(b)
    if fa * fb > 0
        error("f[a] and f[b] should have different signs!")
    end

    c = a
    fc = fa
    c = a
    d = b - a
    e = d

    iter = 0
    maxiter = 10000

    while iter < maxiter
        iter += 1

        if fb * fc > 0
            c = a
            fc = fa
            d = b - a
            e = d
        end

        if abs(fc) < abs(fb)
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        tol = 2.0 * eps() * abs(b) + tol
        m = (c - b) / 2.0

        if abs(m) > tol && abs(fb) > 0
            if abs(e) < tol || abs(fa) <= abs(fb)
                d = m
                e = m
            else
                s = fb / fa
                if a == c
                    p = 2 * m * s
                    q = 1 - s
                else
                    q = fa / fc
                    r = fb / fc
                    p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                    q = (q - 1) * (r - 1) * (s - 1)
                end
                if p > 0
                    q = -q
                else
                    p = -p
                end
                s = e
                e = d
                if 2 * p < 3 * m * q - abs(tol * q) && (p < abs(s * q / 2))
                    d = p / q
                else
                    d = m
                    e = m
                end
            end
            a = b
            fa = fb
            if abs(d) > tol
                b = b + d
            else
                if m > 0
                    b = b + tol
                else
                    b = b - tol
                end
            end
        else
            break
        end
        fb = f(b)
    end
    return b, iter
end

@doc raw"""
    CustomBrent(f::Function, a::Real, b::Real; tol = 1e-14)

Find the root of a function using a customized version of Brent's method. This
implementation adapts Brent's method by incorporating initial guesses for value functions
and distributions, based on linear interpolations from the previous two iterations.

The function is designed for cases where the root-finding process needs to account for
additional parameters beyond the function's value at the endpoints.

# Arguments
- `f::Function`: The function for which the root is being found. It should return a tuple,
  where the first element (`f(z)[1]`) represents the function value at `z`, and the
  subsequent elements represent additional parameters used for interpolation.
- `a::Real`: The lower bound of the interval where the root is located.
- `b::Real`: The upper bound of the interval where the root is located.
- `tol::Real`: The tolerance for convergence. The default is `1e-14`.

# Returns
- `b::Real`: The estimated root of the function.
- `iter::Int`: The number of iterations used to find the root.
- `fb::Tuple`: The function value at the root, along with any interpolated parameters.

# Errors
- Throws an error if the function values at `a` and `b` have the same sign, indicating no
  root exists in the specified interval.
"""
function CustomBrent(f::Function, a::Real, b::Real; tol = 1e-14)
    # Implementation of Brent's method to find a root of a function (as on wikipedia)
    fa = f(a)
    fb = f(b, true, fa[2], fa[3], fa[4])
    if fa[1] * fb[1] > 0
        error("f[a] and f[b] should have different signs!")
    end

    c = a
    fc = fa   # at the beginning: c = a
    c = a
    d = b - a
    e = d

    iter = 0
    maxiter = 10000
    initial = false
    while iter < maxiter
        iter += 1
        if fb[1] * fc[1] > 0
            c = a
            fc = fa
            d = b - a
            e = d
        end

        if abs(fc[1]) < abs(fb[1])
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        tol = 2.0 * eps() * abs(b) + tol
        m = (c - b) / 2.0

        if abs(m) > tol && abs(fb[1]) > 0
            if abs(e) < tol || abs(fa[1]) <= abs(fb[1])
                d = m
                e = m
            else
                s = fb[1] / fa[1]
                if a == c
                    p = 2 * m * s
                    q = 1 - s
                else
                    q = fa[1] / fc[1]
                    r = fb[1] / fc[1]
                    p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                    q = (q - 1) * (r - 1) * (s - 1)
                end
                if p > 0
                    q = -q
                else
                    p = -p
                end
                s = e
                e = d
                if 2 * p < 3 * m * q - abs(tol * q) && (p < abs(s * q / 2))
                    d = p / q
                else
                    d = m
                    e = m
                end
            end
            a = b
            fa = fb
            if abs(d) > tol
                b = b + d
            else
                if m > 0
                    b = b + tol
                else
                    b = b - tol
                end
            end
        else
            break
        end
        fb = f(
            b,
            initial,
            (fa[2] .+ d / (c - a) .* (fc[2] - fa[2])),
            (fa[3] .+ d / (c - a) .* (fc[3] - fa[3])),
            (fa[4] .+ d / (c - a) .* (fc[4] - fa[4])),
        )
    end
    return b, iter, fb
end
