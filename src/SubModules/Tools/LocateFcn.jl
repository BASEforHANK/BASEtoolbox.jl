"""
    locate(x::Number, xx::AbstractVector)

Wrapper around `exp_search`, see its documentation for details.
"""
locate(x::Number, xx::AbstractVector) = exp_search(x, xx)

"""
    bin_search(x::Number, xx::AbstractVector)

Perform a binary search to locate the position of `x` in the sorted vector `xx`.

# Arguments

  - `x::Number`: The value whose nearest neighbors are to be found.
  - `xx::AbstractVector`: A sorted vector serving as the lookup table.

# Returns

  - `j::Int`: The index such that `xx[j] ≤ x < xx[j+1]`. Returns `1` if `x` is less than or
    equal to the first element of `xx`, `length(xx)` if `x` is greater than or equal to the
    last element of `xx`.

# Notes

  - The input vector `xx` must be sorted in ascending order.
"""
function bin_search(x::Number, xx::AbstractVector)
    N = length(xx)

    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        jl = 1
        ju = N
        while (ju - jl) != 1
            jm = div((ju + jl), 2)
            @inbounds if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end

    return j
end

"""
    exp_search(x::Number, xx::AbstractVector)

Perform an exponential search to locate the position of `x` in the sorted vector `xx`.
Exponential search quickly narrows down the search interval in a first step, followed by a
binary search for precise localization in a second step.

# Arguments

  - `x::Number`: The value whose nearest neighbors are to be found.
  - `xx::AbstractVector`: A sorted vector serving as the lookup table.

# Returns

  - `j::Int`: The index such that `xx[j] ≤ x < xx[j+1]`. Returns `1` if `x` is less than or
    equal to the first element of `xx`, `length(xx)` if `x` is greater than or equal to the
    last element of `xx`.

# Notes

  - The input vector `xx` must be sorted in ascending order.
"""
function exp_search(x::Number, xx::AbstractVector)
    N = length(xx)
    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        bound = 2
        @inbounds while bound < N && x > xx[bound]
            bound *= 2
        end
        jl = div(bound, 2)
        ju = min(N, bound)
        while (ju - jl) != 1
            jm = div((ju + jl), 2)
            @inbounds if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end

    return j
end
