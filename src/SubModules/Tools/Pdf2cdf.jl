"""
    cdf_to_pdf(cdf)

Computes the probability density function (PDF) from a given cumulative distribution
function (CDF) by taking finite differences along all dimensions.

# Arguments

  - `cdf`: An array representing the cumulative distribution function. Can be
    one-dimensional or higher-dimensional.

# Returns

  - `Y`: An array of the same size as `cdf`, representing the computed probability density
    function.

# Notes

Allows for one-asset and two-asset case by only taking differences where there are at least
two elements.
"""
function cdf_to_pdf(cdf)
    Y = copy(cdf)
    for j = 1:length(size(cdf))
        if size(cdf, j) > 1
            t = collect(size(cdf))
            t[j] = 1
            aux = zeros(eltype(Y), Tuple(t))
            Y = diff(cat([aux, Y]...; dims = j); dims = j)
        end
    end
    return Y
end

"""
    pdf_to_cdf(pdf)

Computes the cumulative distribution function (CDF) from a given probability density
function (PDF) by taking cumulative sums along all dimensions.

# Arguments

  - `pdf`: An array representing the probability density function. Can be one-dimensional or
    higher-dimensional.

# Returns

  - `Y`: An array of the same size as `pdf`, representing the computed cumulative
    distribution function.
"""
function pdf_to_cdf(pdf)
    Y = copy(pdf)
    for j = 1:length(size(pdf))
        Y = cumsum(Y; dims = j)
    end
    return Y
end
