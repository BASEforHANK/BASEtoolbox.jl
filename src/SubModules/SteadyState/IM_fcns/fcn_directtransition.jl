"""
    DirectTransition(
        b_a_star::Array,
        b_n_star::Array,
        k_a_star::Array,
        distr::Array,
        λ,
        Π::Array,
        n_par,
    )

Function to calculate next periods distribution of households given todays policy functions
and the current distribution.

The function calls the in-place version [`DirectTransition!`](@ref) to update the
distribution. For details on the implementation see that function and subsection
'Aggregation via non-stochastic simulations in ['Computational Notes.md'](Computational
Notes.md) for further details.

# Arguments

  - `b_a_star::Array`: Liquid asset policy function for adjustment case
  - `b_n_star::Array`: Liquid asset policy function for non-adjustment case
  - `k_a_star::Array`: Illiquid asset policy function for adjustment case
  - `distr::Array`: Current distribution of households
  - `λ::Float64`: Probability of adjustment
  - `Π::Array`: Transition matrix for idiosyncratic productivity states
  - `n_par::NumericalParameters`: Struct holding numerical parameters

# Returns

  - `dPrime`: Updated distribution of households
"""
function DirectTransition(
    b_a_star::Array,
    b_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par,
)
    dPrime = zeros(eltype(distr), size(distr))
    DirectTransition!(dPrime, b_a_star, b_n_star, k_a_star, distr, λ, Π, n_par, n_par.model)
    return dPrime
end

"""
    DirectTransition!(
        dPrime,
        b_a_star::Array,
        b_n_star::Array,
        k_a_star::Array,
        distr::Array,
        λ,
        Π::Array,
        n_par,
        model::Union{OneAsset, TwoAsset, CompleteMarkets},
    )

Function to calculate next periods distribution of households given todays policy functions
and the current distribution. The function updates the distribution in place.

The function is used in `FSYS.jl` function to update the distribution of households given
the policy functions and the current distribution. For details on the implementation see
subsection 'Aggregation via non-stochastic simulations in ['Computational
Notes.md'](Computational Notes.md) for further details.

# Arguments

  - `dPrime::Array`: To be updated distribution of households
  - `b_a_star::Array`: Liquid asset policy function for adjustment case
  - `b_n_star::Array`: Liquid asset policy function for non-adjustment case
  - `k_a_star::Array`: Illiquid asset policy function for adjustment case
  - `distr::Array`: Current distribution of households
  - `λ::Float64`: Probability of adjustment
  - `Π::Array`: Transition matrix for idiosyncratic productivity states
  - `n_par::NumericalParameters`: Struct holding numerical parameters
  - `model`: Model type, either `CompleteMarkets`, `OneAsset`, or `TwoAsset`

# Returns

  - `dPrime`: Updated distribution of households
"""
function DirectTransition!(
    dPrime,
    b_a_star::Array,
    b_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par,
    model::TwoAsset,
)
    # Create linear interpolation weights from policy functions for both cases
    j_k_a, ω_right_k_a = MakeWeightsLight(k_a_star, n_par.grid_k)
    j_b_a, ω_right_b_a = MakeWeightsLight(b_a_star, n_par.grid_b)
    j_b_n, ω_right_b_n = MakeWeightsLight(b_n_star, n_par.grid_b)

    # Setup the blockindex necessary for linear indexing over productivity state tomorrow
    # Below we will use this to calculate the target index of the policy functions for the
    # transition matrix. We need to add the blockindex to the linear index of the policy
    # function to get the correct target index associated with todays choices as function of
    # the state space and tomorrows productivity state
    blockindex = (0:(n_par.nh - 1)) * n_par.nk * n_par.nb

    # Loop over all states of the current period and calculate the updated distribution
    @inbounds begin
        for hh = 1:(n_par.nh) # all current productivity states
            for kk = 1:(n_par.nk) # all current illiquid asset states
                for bb = 1:(n_par.nb) # all current liquid asset states

                    # current mass of households in this state
                    dd = distr[bb, kk, hh]

                    # linear index of the policy functions for the current state
                    j_adj = (j_b_a[bb, kk, hh] .+ (j_k_a[bb, kk, hh] .- 1) .* n_par.nb)
                    j_non = (j_b_n[bb, kk, hh] .+ (kk - 1) .* n_par.nb)

                    # future mass for the liquid assets of non adjusters
                    ω = ω_right_b_n[bb, kk, hh]
                    d_L_n = (1.0 .- λ) .* (dd .* (1.0 .- ω))
                    d_R_n = (1.0 .- λ) .* (dd .* ω)

                    # future mass for the illiquid assets of adjusters
                    ω = ω_right_k_a[bb, kk, hh]
                    d_L_k_a = λ .* (dd .* (1.0 .- ω))
                    d_R_k_a = λ .* (dd .* ω)

                    # future mass for both assets of adjusters
                    ω = ω_right_b_a[bb, kk, hh]
                    d_LL_a = (d_L_k_a .* (1.0 .- ω))
                    d_LR_a = (d_L_k_a .* ω)
                    d_RL_a = (d_R_k_a .* (1.0 .- ω))
                    d_RR_a = (d_R_k_a .* ω)

                    # transitions to future productivity states
                    for hh_prime = 1:(n_par.nh)
                        # extract the probability of transitioning from the current
                        # productivity state to the future productivity state
                        pp = Π[hh, hh_prime]

                        # linear index of the policy functions for the future state
                        j_a = j_adj .+ blockindex[hh_prime]
                        j_n = j_non .+ blockindex[hh_prime]

                        # update the distribution of households by adding the weighted mass
                        # of households that will transition to the future state
                        # corresponding to the gammas in the transition matrix
                        dPrime[j_a] += pp .* d_LL_a
                        dPrime[j_a + 1] += pp .* d_LR_a
                        dPrime[j_a + n_par.nb] += pp .* d_RL_a
                        dPrime[j_a + n_par.nb + 1] += pp .* d_RR_a
                        dPrime[j_n] += pp .* d_L_n
                        dPrime[j_n + 1] += pp .* d_R_n
                    end
                end
            end
        end
    end
end

function DirectTransition!(
    dPrime,
    b_a_star::Array,
    b_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par,
    model::OneAsset,
)
    # Create linear interpolation weights from policy functions for both cases
    j_b_n, ω_right_b_n = MakeWeightsLight(b_n_star, n_par.grid_b)

    # Setup the blockindex necessary for linear indexing over productivity state tomorrow
    # Below we will use this to calculate the target index of the policy functions for the
    # transition matrix. We need to add the blockindex to the linear index of the policy
    # function to get the correct target index associated with todays choices as function of
    # the state space and tomorrows productivity state
    blockindex = (0:(n_par.nh - 1)) * n_par.nk * n_par.nb

    # Loop over all states of the current period and calculate the updated distribution
    @inbounds begin
        for hh = 1:(n_par.nh) # all current productivity states
            for kk = 1:(n_par.nk) # all current illiquid asset states
                for bb = 1:(n_par.nb) # all current liquid asset states

                    # current mass of households in this state
                    dd = distr[bb, kk, hh]

                    # linear index of the policy functions for the current state
                    j_non = (j_b_n[bb, kk, hh] .+ (kk - 1) .* n_par.nb)

                    # future mass for the liquid assets of non adjusters
                    ω = ω_right_b_n[bb, kk, hh]
                    d_L_n = (dd .* (1.0 .- ω))
                    d_R_n = (dd .* ω)

                    # transitions to future productivity states
                    for hh_prime = 1:(n_par.nh)
                        # extract the probability of transitioning from the current
                        # productivity state to the future productivity state
                        pp = Π[hh, hh_prime]

                        # linear index of the policy functions for the future state
                        j_n = j_non .+ blockindex[hh_prime]

                        # update the distribution of households by adding the weighted mass
                        # of households that will transition to the future state
                        # corresponding to the gammas in the transition matrix
                        dPrime[j_n] += pp .* d_L_n
                        dPrime[j_n + 1] += pp .* d_R_n
                    end
                end
            end
        end
    end
end

function DirectTransition!(
    dPrime,
    b_a_star::Array,
    b_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par,
    model::CompleteMarkets,
)
    # don't do anything
end
