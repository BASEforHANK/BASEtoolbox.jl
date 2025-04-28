"""
    centraldiff(x)

Calculate the central difference of a vector `x`. The central difference is a numerical
approximation to the derivative of a function, computed using the values of the function at
adjacent points. The first and last elements of the vector are computed using forward and
backward differences, respectively, because a central difference cannot be computed at the
boundaries.

# Arguments

  - `x::AbstractVector`: The input vector for which the central difference is computed.

# Returns

  - `a::AbstractVector`: The central difference approximation of the input vector `x`. The
    output vector has the same length as `x`, with the first and last elements computed
    using forward and backward differences, respectively.

# Errors:

  - Throws an error if the input x is not a 1-dimensional vector.
"""
function centraldiff(x)
    if length(x) == 1
        return 0.0
    end
    if ndims(x) != 1
        error("wrong number of dimensions, only 1 allowed")
    else
        a = similar(x)
        dx = diff(x) ./ 2.0
        a .= [dx; dx[end]]
        a .+= [dx[1]; dx]
        return a
    end
end

"""
    centralderiv(y, x, dims)

Calculate the approximate derivative of the array `y` with respect to the array `x` along
the specified dimension(s) `dims` using central differences (cf.
[`Tools.centraldiff()`](@ref)).

# Arguments

  - `y::AbstractArray`: The array for which the derivative is computed.
  - `x::AbstractArray`: The array with respect to which the derivative of `y` is calculated.
  - `dims::Union{Int, Tuple{Vararg{Int}}}`: The dimension(s) along which to calculate the
    derivative.

# Returns

  - `AbstractArray`: The approximate derivative of `y` with respect to `x`.
"""
function centralderiv(y, x, dims)
    dy = mapslices(centraldiff, y; dims = dims)
    dx = mapslices(centraldiff, x; dims = dims) .+ eps()
    return dy ./ dx
end
