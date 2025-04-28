@doc raw"""
    Fastroot(xgrid::Vector, fx::AbstractArray)

Used for fast root finding using linear interpolation and a simplified Newton step.

This function performs a fast root-finding algorithm on a given grid (`xgrid`) and function
values (`fx`). It applies a single Newton step at the largest negative function value.

# Arguments
- `xgrid::Vector`: The grid of x values.
- `fx::AbstractArray`: Function values evaluated at the points in `xgrid`. The array should
  be reshaped to match the dimensions of `xgrid`.

# Returns
- `roots::Vector`: A vector of roots found through interpolation for each column of `fx`.
  The length of the vector corresponds to the number of columns in `fx`.
"""
function Fastroot(xgrid::Vector, fx::AbstractArray)
    #fast linear interpolation root finding (=one Newton step at largest negative function
    #value) stripped down version of interp1 that accepts multiple inputs [max 3] that are
    #   interpolated over the same grids x & xi
    xgrid = xgrid[:]
    fx = reshape(fx, length(xgrid), :)
    # Make use of the fact that the difference equation is monotonically increasing in m,
    # use sign for robustness.
    n = size(fx)
    roots = Array{eltype(fx),1}(undef, n[2])
    fz = similar(fx[:, 1])
    sfz = sign.(fz)

    @views @inbounds begin
        for j = 1:n[2]
            fz = fx[:, j]
            sfz .= sign.(fz)
            idx = locate(0.0, sfz)
            if idx >= n[1]
                roots[j] = xgrid[end]
            elseif idx <= 1
                roots[j] = xgrid[1]
            else
                fxx = fz[idx]
                dfxx = fz[idx .+ 1] .- fxx
                xx = xgrid[idx]
                dxx = xgrid[idx .+ 1] .- xgrid[idx]
                roots[j] = xx .- fxx .* dxx ./ dfxx
            end
        end
    end
    return roots
end
