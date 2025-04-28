"""
    MakeWeights(xpol, grid)

This function takes a vector of points `xpol` and a grid `grid` and returns the index of the
grid point to the left of each point in `xpol`, the weight associated with this grid point
to the left of each point in `xpol`, and the weight associated with the grid point to the
right of each point in `xpol`.

The function is used in the [`MakeTransition()`](@ref) and
[`MultipleDirectTransition`](@ref) functions to compute the weights for the gambles implied
by the policy functions when applying the method of [Young
(2010)](https://www.sciencedirect.com/science/article/abs/pii/S0165188909001316). Refer to
subsection 'Aggregation via non-stochastic simulations in ['Computational
Notes.md'](Computational Notes.md) for further details.

# Arguments

  - `xpol::Array{Float64, 1}`: Vector of points to be interpolated
  - `grid::Array{Float64, 1}`: Grid to interpolate on

# Returns

  - `idx::Array{Int64, 1}`: Index of the grid point to the left of each point in `xpol`
  - `weightleft::Array{Float64, 1}`: Weight associated with the grid-point to the left of
    each point in `xpol`, the higher the closer the point is to the grid-point with index
    `idx`
  - `weightright::Array{Float64, 1}`: Weight associated with the grid-point to the right of
    each point in  `xpol`

# Example

```julia
xpol = [0.1, 0.2, 0.3, 0.4, 0.5]
grid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
idx, weightleft, weightright = MakeWeights(xpol, grid)
```
"""
function MakeWeights(xpol, grid)

    # containers for outputs
    idx = zeros(Int64, size(xpol))
    weightright = zeros(typeof(xpol[1]), size(xpol))
    weightleft = zeros(typeof(xpol[1]), size(xpol))

    # calculate the distance between grid points
    dx = diff(grid)

    # iterate over the entries in the policy function
    for i in eachindex(xpol)
        if xpol[i] .<= grid[1]
            # check if the point is below the lowest grid point set the index to the first
            # grid point effectively allowing for extrapolation below the grid
            idx[i] = 1

        elseif xpol[i] .>= grid[end]
            # check if the point is above the highest grid point set the index to the last
            # grid point effectively allowing for extrapolation above the grid
            idx[i] = length(grid) .- 1

        else
            # if the point is internal find the index of the grid point to the left of the
            # point
            idx[i] = locate(xpol[i], grid)
        end

        # calculate the associated weight to the grid point to the right and left of the
        # point according to the formula for ω_k and ω_b in the computational notes and
        # following Young (2010)
        weightright[i] = (xpol[i] - grid[idx[i]]) / dx[idx[i]]
        weightleft[i] = 1.0 - weightright[i]
    end
    return idx, weightleft, weightright
end

"""
    MakeWeightsLight(xpol, grid)

This function takes a vector of points `xpol` and a grid `grid` and returns the index of the
grid point to the left of each point in `xpol` and the weight to the right of each point in
`xpol`. It is an optimized version of the [`MakeWeights`](@ref) function.

The function is used in the [`DirectTransition`](@ref) function to calculate weights from
policy functions to compute the weights for the gambles when applying the method of [Young
(2010)](https://www.sciencedirect.com/science/article/abs/pii/S0165188909001316). Refer to
subsection 'Aggregation via non-stochastic simulations in ['Computational
Notes.md'](Computational Notes.md) for further details.

# Arguments

  - `xpol::Array{Float64, 1}`: Vector of points to be interpolated
  - `grid::Array{Float64, 1}`: Grid to interpolate on

# Returns

  - `idx::Array{Int64, 1}`: Index of the grid point to the left of each point in `xpol`
  - `weightright::Array{Float64, 1}`: Weight to the right of each point in `xpol`

# Example

```julia
xpol = [0.1, 0.2, 0.3, 0.4, 0.5]
grid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
idx, weightright = MakeWeightsLight(xpol, grid)
```
"""
function MakeWeightsLight(xpol, grid)

    # containers for outputs
    idx = zeros(Int64, size(xpol))
    weightright = zeros(eltype(xpol), size(xpol))

    # calculate the distance between grid points
    dx = diff(grid)

    # iterate over the entries in the policy function
    @inbounds begin
        for i in eachindex(xpol)
            if xpol[i] .<= grid[1]
                # check if the point is below the lowest grid point set the index to the
                # first grid point effectively allowing for extrapolation below the grid
                idx[i] = 1

            elseif xpol[i] .>= grid[end]
                # check if the point is above the highest grid point set the index to the
                # last grid point effectively allowing for extrapolation above the grid
                idx[i] = length(grid) .- 1

            else
                # the point is internal find the index of the grid point to the left of the
                # point
                idx[i] = locate(xpol[i], grid)
            end

            # calculate the weight associated with the grid point to the right of the point
            # according to the formula for (1-ω_k) and (1-ω_b) in the computational notes
            # and following Young (2010) The higher weightright the closer the grid-point to
            # the right of the point is to the point.
            weightright[i] = (xpol[i] .- grid[idx[i]]) ./ dx[idx[i]]
        end
    end
    return idx, weightright
end
