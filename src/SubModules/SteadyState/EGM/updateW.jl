"""
    updateW(
        EWkPrime::Array,
        x_a_star::Array,
        x_n_star::Array,
        b_n_star::Array,
        args_hh_prob::Vector,
        m_par,
        n_par,
        model::Union{CompleteMarkets,OneAsset,TwoAsset},
    )

Update the marginal values for the liquid and illiquid assets, given the current
continuation value for the illiquid asset `EWkPrime`, the policy functions [`x_a_star`,
`x_n_star`, `b_n_star`], and today's prices [`args_hh_prob`].

Consult the document ['Computational Notes.md'](Computational Notes.md), specifically
section 3 with the title 'Update the continuation values (CV1) and (CV2) using the Envelope
conditions' for an explanation of the function's code.

# Arguments

  - `EWkPrime`: Marginal continuation values of illiquid assets
  - `x_a_star`: Optimal (on-grid) policy for the composite with adjustment of illiquid asset
  - `x_n_star`: Optimal (on-grid) policy for the composite without adjustment of illiquid
    asset
  - `b_n_star`: Optimal (on-grid) policy for liquid asset without adjustment of illiquid
    asset
  - `args_hh_prob`: Vector of arguments to the household problem
  - `n_par`, `m_par`
  - `model`: Model type, either `CompleteMarkets`, `OneAsset`, or `TwoAsset`

# Returns

  - `Wk`: Marginal value of illiquid assets
  - `Wb`: Marginal value of liquid assets
"""
function updateW(
    EWkPrime::Array,
    x_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    args_hh_prob::Vector,
    m_par,
    n_par,
    model::Union{CompleteMarkets,OneAsset,TwoAsset},
)

    # Preallocate variables
    Wk = similar(EWkPrime)
    Wb = similar(EWkPrime)
    mutil_x_a = similar(EWkPrime)

    # Call the inplace function to update the marginal value functions
    updateW!(
        Wk,
        Wb,
        mutil_x_a,
        EWkPrime,
        x_a_star,
        x_n_star,
        b_n_star,
        args_hh_prob,
        m_par,
        n_par,
        model,
    )

    # Return results
    return Wk, Wb
end

"""
    updateW!(
        Wk::Array,
        Wb::Array,
        mutil_x_a::Array,
        EWkPrime::Array,
        x_a_star::Array,
        x_n_star::Array,
        b_n_star::Array,
        args_hh_prob::Vector,
        m_par,
        n_par,
        model::Union{CompleteMarkets,OneAsset,TwoAsset},
    )

In-place version of [`updateW()`](@ref), see that function for details.
"""
function updateW!(
    Wk::Array,
    Wb::Array,
    mutil_x_a::Array,
    EWkPrime::Array,
    x_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    args_hh_prob::Vector,
    m_par,
    n_par,
    model::CompleteMarkets,
)
    # Does nothing, as the model is complete markets
    Wb = ones(n_par.nb, n_par.nk, n_par.nh)
    Wk = ones(n_par.nb, n_par.nk, n_par.nh)
end

function updateW!(
    Wk::Array,
    Wb::Array,
    mutil_x_a::Array,
    EWkPrime::Array,
    x_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    args_hh_prob::Vector,
    m_par,
    n_par,
    model::OneAsset,
)
    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Unpack
    ## ------------------------------------------------------------------------------------

    nb, nk, nh = size(x_n_star)
    β::Float64 = m_par.β

    ## ------------------------------------------------------------------------------------
    ## Update marginal value of liquid assets using the Envelope conditions
    ## ------------------------------------------------------------------------------------

    #=
    Compute marginal utility of the composite, non-adjustment case.

    Note: Here we use Wb as a temporary storage for the marginal utility of the composite.
    Later we overwrite Wb with the probability-weighted sum of the adjustment and
    non-adjustment cases of the marginal value of liquid assets.

    For the probability weihgted sum of the marginal value of liquid assets, see the
    documentation of algorithm (CV2) with plugged in Envelope Conditions (EC2) and (EC3).

    Note: The marginal value of liquid assets is adjusted for the VAT rate but not yet
    premultiplied by real interest rates for technical reasons below. For EWb to fulfill
    (CV2), effective interest rates are multiplied outside of this function.
    =#

    # Compute marginal utility of the composite, non-adjustment case
    mutil!(Wb, x_n_star, m_par)

    Wk .= 1.0
    mutil_x_a .= 0.0
end

function updateW!(
    Wk::Array,
    Wb::Array,
    mutil_x_a::Array,
    EWkPrime::Array,
    x_a_star::Array,
    x_n_star::Array,
    b_n_star::Array,
    args_hh_prob::Vector,
    m_par,
    n_par,
    model::TwoAsset,
)
    @read_args_hh_prob()

    ## ------------------------------------------------------------------------------------
    ## Unpack
    ## ------------------------------------------------------------------------------------

    nb, nk, nh = size(x_n_star)
    β::Float64 = m_par.β

    ## ------------------------------------------------------------------------------------
    ## Update marginal value of liquid assets using the Envelope conditions
    ## ------------------------------------------------------------------------------------

    #=
    Compute marginal utility of the composite, non-adjustment case.

    Note: Here we use Wb as a temporary storage for the marginal utility of the composite.
    Later we overwrite Wb with the probability-weighted sum of the adjustment and
    non-adjustment cases of the marginal value of liquid assets.

    For the probability weihgted sum of the marginal value of liquid assets, see the
    documentation of algorithm (CV2) with plugged in Envelope Conditions (EC2) and (EC3).

    Note: The marginal value of liquid assets is adjusted for the VAT rate but not yet
    premultiplied by real interest rates for technical reasons below. For EWb to fulfill
    (CV2), effective interest rates are multiplied outside of this function.
    =#

    # Compute marginal utility of the composite, non-adjustment case
    mutil!(Wb, x_n_star, m_par)

    # Compute marginal utility of the composite, adjustment case
    mutil!(mutil_x_a, x_a_star, m_par)

    # Probability-weighted sum of marginal value of liquid assets
    Wb .= m_par.λ .* mutil_x_a .+ (1.0 - m_par.λ) .* Wb
    Wb .= Wb ./ (1.0 .+ (Tc .- 1.0))

    ## ------------------------------------------------------------------------------------
    ## Update marginal value of illiquid assets
    ## ------------------------------------------------------------------------------------

    #=
    Interpolation of marginal value of illiquid assets over the savings function of
    non-adjustment case; inverse here to have a strictly increasing function to allow for
    extrapolatation of a strictly decreasing function EWkPrime (to avoid negative values)

    The interpolated values are allocated to Wk, which is later overwritten by the
    probability-weighted sum of the adjustment and non-adjustment cases
    =#
    @inbounds @views begin
        for j::Int = 1:nh
            for k::Int = 1:nk
                mylinearinterpolate!(
                    Wk[:, k, j],
                    n_par.grid_b,
                    invmutil(EWkPrime[:, k, j], m_par),
                    b_n_star[:, k, j],
                )
            end
        end
    end

    # Revert the inversion
    Wk .= mutil(Wk, m_par)

    # Probability-weighted sum of marginal value of illiquid assets in both cases of
    # adjustment, see documentation of algorithm (CV1a) for details
    Wk .=
        (RK .- 1.0) .* Wb .+ m_par.λ .* q .* mutil_x_a ./ (1.0 .+ (Tc .- 1.0)) .+
        (1.0 .- m_par.λ) .* β .* Wk
end
