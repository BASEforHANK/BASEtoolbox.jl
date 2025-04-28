"""
    tot_dual(x::ForwardDiff.Dual)

Computes the total derivative of a dual number by summing over all its partial derivatives.

# Arguments

  - `x::ForwardDiff.Dual`: A dual number from `ForwardDiff.jl`.

# Returns

  - `a`: The sum of all partial derivatives of `x`.
"""
function tot_dual(x::ForwardDiff.Dual)
    a = sum(ForwardDiff.partials(x, :))
    return a
end

"""
    realpart(x::ForwardDiff.Dual)

Extracts the real (primal) part of a dual number. Multiple dispatch allows using one
function, `realpart()`, with both types, `ForwardDiff.Dual` and `Float64`.

# Arguments

  - `x::ForwardDiff.Dual`: A dual number from `ForwardDiff.jl`.

# Returns

  - `a`: The real (primal) part of the dual number `x`, extracted using `ForwardDiff.value`.
"""
function realpart(x::ForwardDiff.Dual)
    a = ForwardDiff.value(x)
    return a
end

"""
    realpart(x::Float64)

Returns the input value `x`. Multiple dispatch allows using one function, `realpart()`, with
both types, `ForwardDiff.Dual` and `Float64`.

# Arguments

  - `x::Float64`: A real number of type `Float64`.

# Returns

  - `a`: The input value `x`.
"""
function realpart(x::Float64)
    a = x
    return a
end

"""
    dualpart(x::ForwardDiff.Dual)

Extracts the dual (derivative) part of a dual number. Multiple dispatch allows using one
function, `dualpart()`, with both types, `ForwardDiff.Dual` and `Float64`.

# Arguments

  - `x::ForwardDiff.Dual`: A dual number from `ForwardDiff.jl`.

# Returns

  - `b`: The dual (derivative) part of the dual number `x`, extracted from the partial
    derivatives.
"""
function dualpart(x::ForwardDiff.Dual)
    a = ForwardDiff.partials.(x)
    b = a.values[1]

    return b
end

"""
    dualpart(x::Float64)

Returns the input value `x`. Multiple dispatch allows using one function, `dualpart()`, with
both types, `ForwardDiff.Dual` and `Float64`.

# Arguments

  - `x::Float64`: A real number of type `Float64`.

# Returns

  - `a`: The input value `x`.
"""
function dualpart(x::Float64)
    b = x

    return b
end
