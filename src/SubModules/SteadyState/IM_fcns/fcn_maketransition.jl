"""
    function MakeTransition(
        b_a_star::Array{Float64,3},
        b_n_star::Array{Float64,3},
        k_a_star::Array{Float64,3},
        Π::Array{Float64,2},
        n_par,
        model::Union{OneAsset,TwoAsset},
    )

This function is used to calculate the inputs to generate the sparse stochastic transition
matrix Γ of households over the state space. The transition function is calculated in
[`fcn_ksupply()`](@ref).

Calculate the mass (weight) of households that transition from the startindex to the
targetindex according to the method of [Young
(2010)](https://www.sciencedirect.com/science/article/abs/pii/S0165188909001316). Refer to
subsection 'Aggregation via non-stochastic simulations' in ['Computational
Notes.md'](Computational Notes.md) for further details.

Requires global functions from the [`Tools`](@ref) submodule.

# Arguments

  - `b_a_star::Array{Float64, 3}`: Liquid asset policy function for adjustment case
  - `b_n_star::Array{Float64, 3}`: Liquid asset policy function for non-adjustment case
  - `k_a_star::Array{Float64, 3}`: Illiquid asset policy function for adjustment case
  - `Π::Array{Float64, 2}`: Transition matrix for idiosyncratic productivity states
  - `n_par::NumericalParameters`: Struct holding numerical parameters
  - `model`: Model type, either `OneAsset`, or `TwoAsset`

# Returns

  - `S_a`, `S_n`: matrices of integers indicating starting location for creation of `sparse`
    transition matrix Γ (with [`a`] or without [`n`] adjustment of illiquid asset)
  - `T_a`, `T_n`: matrices of intergers indicating target location for creation of `sparse`
    transition matrix Γ (with [`a`] or without [`n`] adjustment of illiquid asset)
  - `W_a`, `W_n`: matrices of weights associated with the probability that a household
    transitions from the starting location to the target location (with [`a`] or without
    [`n`] adjustment of illiquid asset)
"""
function MakeTransition(
    b_a_star::Array{Float64,3},
    b_n_star::Array{Float64,3},
    k_a_star::Array{Float64,3},
    Π::Array{Float64,2},
    n_par,
    model::TwoAsset,
)

    # Create linear interpolation weights from policy functions for both cases and find the
    # associated indices j next on the left of the policy function
    j_k_a, ω_left_k_a, ω_right_k_a = MakeWeights(k_a_star, n_par.grid_k)
    j_b_a, ω_left_b_a, ω_right_b_a = MakeWeights(b_a_star, n_par.grid_b)
    j_b_n, ω_left_b_n, ω_right_b_n = MakeWeights(b_n_star, n_par.grid_b)

    # Create the containers of the output for the adjustment case
    weight_adj =
        Array{eltype(k_a_star),3}(undef, 4, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    targetindex_adj = Array{Int,3}(undef, 4, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    startindex_adj = Array{Int,3}(undef, 4, n_par.nh, n_par.nk * n_par.nb * n_par.nh)

    # Set up the blockindex necessary for linear indexing over productivity state tomorrow
    # Below we will use this to calculate the target index of the policy functions for the
    # transition matrix. We need to add the blockindex to the linear index of the policy
    # function to get the correct target index associated with todays choices as function of
    # the state space and tomorrows productivity state
    blockindex = (0:(n_par.nh - 1)) * n_par.nk * n_par.nb

    # Set the runindex, which is used as startindex_adj
    runindex = 0

    # Loop over all states and calculate weights, targetindex_adj and startindex_adj
    for hh = 1:(n_par.nh) # all current income states
        for kk = 1:(n_par.nk) # all current illiquid asset states
            for bb = 1:(n_par.nb) # all current liquid asset states
                # update the run index by one iteration
                runindex = runindex + 1

                # Calculate the weights for the four possible transitions ωLL: left liquid,
                # left illiquid
                ωLL = ω_left_b_a[bb, kk, hh] .* ω_left_k_a[bb, kk, hh]

                # ωRL: right liquid, left illiquid
                ωRL = ω_right_b_a[bb, kk, hh] .* ω_left_k_a[bb, kk, hh]

                # ωLR: left liquid, right illiquid
                ωLR = ω_left_b_a[bb, kk, hh] .* ω_right_k_a[bb, kk, hh]

                # ωRR: right liquid, right illiquid
                ωRR = ω_right_b_a[bb, kk, hh] .* ω_right_k_a[bb, kk, hh]

                # Calculate the target index, which policy functions point towards, using
                # linear indexing Note: The index is calculated by adding the linear index
                # of the liquid asset policy function to the linear index of the illiquid
                # asset policy function times the number of liquid asset states
                j_adj = j_b_a[bb, kk, hh] .+ (j_k_a[bb, kk, hh] - 1) .* n_par.nb

                # iterating over all future income states
                for jj = 1:(n_par.nh)
                    # extract the probability of transitioning from the current income state
                    # to the future income state
                    pp = Π[hh, jj]

                    # extracting the block index associated with the productivity state
                    bb = blockindex[jj]

                    # calculate the weights associated with all four asset transitions
                    weight_adj[1, jj, runindex] = ωLL .* pp
                    weight_adj[2, jj, runindex] = ωRL .* pp
                    weight_adj[3, jj, runindex] = ωLR .* pp
                    weight_adj[4, jj, runindex] = ωRR .* pp

                    # calculate the target index associated with all four asset transitions
                    # Note: The target index needs to match the linear index of the policy
                    # function that both policy functions are on the left grid point (first
                    # case), both policy funciton are on the right grid point (last case),
                    # and that the policy functions are on different grid points (second and
                    # third case)
                    targetindex_adj[1, jj, runindex] = j_adj .+ bb
                    targetindex_adj[2, jj, runindex] = j_adj + 1 .+ bb
                    targetindex_adj[3, jj, runindex] = j_adj + n_par.nb .+ bb
                    targetindex_adj[4, jj, runindex] = j_adj + n_par.nb + 1 .+ bb

                    # set the start index to the run index
                    startindex_adj[1, jj, runindex] = runindex
                    startindex_adj[2, jj, runindex] = runindex
                    startindex_adj[3, jj, runindex] = runindex
                    startindex_adj[4, jj, runindex] = runindex
                end
            end
        end
    end

    # Flatten the arrays
    S_a = startindex_adj[:]
    T_a = targetindex_adj[:]
    W_a = weight_adj[:]

    # Create the containers of the output for the non adjustment case
    weight_non = zeros(eltype(k_a_star), 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    targetindex_non = zeros(Int, 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    startindex_non = zeros(Int, 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)

    # reset the runindex to 0
    runindex = 0

    # Loop over all states and calculate weights, targetindex and startindex
    for hh = 1:(n_par.nh) # all current income states
        for kk = 1:(n_par.nk) # all current illiquid asset states
            for bb = 1:(n_par.nb) # all current liquid asset states
                # update the run index by one iteration
                runindex = runindex + 1

                # Calculate the weights for the two possible transitions
                ωL = ω_left_b_n[bb, kk, hh]
                ωR = ω_right_b_n[bb, kk, hh]

                # Calculate the target index, which policy functions point towards, using
                # linear indexing Note: The index is calculated by adding the linear index
                # of the liquid asset policy function to the linear index of the illiquid
                # asset holding times the number of liquid asset states
                j_non = j_b_n[bb, kk, hh] .+ (kk - 1) .* n_par.nb

                # iterating over all future income states
                for jj = 1:(n_par.nh)
                    # extract the probability of transitioning from the current income state
                    # to the future income state
                    pp = Π[hh, jj]

                    # calculate the weights associated with all two asset transitions
                    weight_non[1, jj, runindex] = ωL .* pp
                    weight_non[2, jj, runindex] = ωR .* pp

                    # calculate the target index associated with all two asset transitions
                    targetindex_non[1, jj, runindex] = j_non .+ blockindex[jj]
                    targetindex_non[2, jj, runindex] = j_non .+ 1 .+ blockindex[jj]

                    # set the start index to the run index
                    startindex_non[1, jj, runindex] = runindex
                    startindex_non[2, jj, runindex] = runindex
                end
            end
        end
    end

    # Flatten the arrays
    S_n = startindex_non[:]
    T_n = targetindex_non[:]
    W_n = weight_non[:]

    return S_a, T_a, W_a, S_n, T_n, W_n
end

function MakeTransition(
    b_a_star::Array{Float64,3},
    b_n_star::Array{Float64,3},
    k_a_star::Array{Float64,3},
    Π::Array{Float64,2},
    n_par,
    model::OneAsset,
)

    # Create linear interpolation weights from policy functions for both cases and find the
    # associated indices j next on the left of the policy function
    j_b_n, ω_left_b_n, ω_right_b_n = MakeWeights(b_n_star, n_par.grid_b)

    # Create the containers of the output for the non adjustment case
    weight_non = zeros(eltype(k_a_star), 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    targetindex_non = zeros(Int, 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)
    startindex_non = zeros(Int, 2, n_par.nh, n_par.nk * n_par.nb * n_par.nh)

    # Set up the blockindex necessary for linear indexing over productivity state tomorrow
    # Below we will use this to calculate the target index of the policy functions for the
    # transition matrix. We need to add the blockindex to the linear index of the policy
    # function to get the correct target index associated with todays choices as function of
    # the state space and tomorrows productivity state
    blockindex = (0:(n_par.nh - 1)) * n_par.nk * n_par.nb

    # reset the runindex to 0
    runindex = 0

    # Loop over all states and calculate weights, targetindex and startindex
    for hh = 1:(n_par.nh) # all current income states
        for kk = 1:(n_par.nk) # all current illiquid asset states
            for bb = 1:(n_par.nb) # all current liquid asset states
                # update the run index by one iteration
                runindex = runindex + 1

                # Calculate the weights for the two possible transitions
                ωL = ω_left_b_n[bb, kk, hh]
                ωR = ω_right_b_n[bb, kk, hh]

                # Calculate the target index, which policy functions point towards, using
                # linear indexing Note: The index is calculated by adding the linear index
                # of the liquid asset policy function to the linear index of the illiquid
                # asset holding times the number of liquid asset states
                j_non = j_b_n[bb, kk, hh] .+ (kk - 1) .* n_par.nb

                # iterating over all future income states
                for jj = 1:(n_par.nh)
                    # extract the probability of transitioning from the current income state
                    # to the future income state
                    pp = Π[hh, jj]

                    # calculate the weights associated with all two asset transitions
                    weight_non[1, jj, runindex] = ωL .* pp
                    weight_non[2, jj, runindex] = ωR .* pp

                    # calculate the target index associated with all two asset transitions
                    targetindex_non[1, jj, runindex] = j_non .+ blockindex[jj]
                    targetindex_non[2, jj, runindex] = j_non .+ 1 .+ blockindex[jj]

                    # set the start index to the run index
                    startindex_non[1, jj, runindex] = runindex
                    startindex_non[2, jj, runindex] = runindex
                end
            end
        end
    end

    # Flatten the arrays
    S_n = startindex_non[:]
    T_n = targetindex_non[:]
    W_n = weight_non[:]

    return S_n, T_n, W_n, S_n, T_n, W_n
end
