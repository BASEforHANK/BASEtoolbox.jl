#=

Template function for Fsys_agg.jl

Explanation:
During the preprocessing step, `PreprocessInputs.jl` reads in this template file and
copies the content from `input_aggregate_model.mod` into the code block below at
the line marked with "# aggregate model marker". The code block is then written to
`FSYS_agg_generated.jl` in the `generated_fcns` directory.

=#

@doc raw"""
    Fsys_agg(X, XPrime, XSS, distrSS, m_par, n_par, indexes)

Return deviations from aggregate equilibrium conditions.

`indexes` can be both `IndexStruct` or `IndexStructAggr`; in the latter case
(which is how function is called by [`LinearSolution_reduced_system()`](@ref)), variable-vectors
`X`,`XPrime`, and `XSS` only contain the aggregate variables of the model.
"""
function Fsys_agg(
    X::AbstractArray,
    XPrime::AbstractArray,
    XSS::Array{Float64,1},
    distrSS::AbstractArray,
    m_par,
    n_par,
    indexes,
)

    ## ------------------------------------------------------------------------------------
    ## Preamble
    ## ------------------------------------------------------------------------------------

    # Initialize the output vector, use the same type as the input
    F = zeros(eltype(X), size(X))

    # Unpack X, XPrime, and XSS into variables
    @generate_equations()

    ## ------------------------------------------------------------------------------------
    ## Get the equilibrium conditions for the aggregate model (`input_aggregate_model.mod`)
    ## ------------------------------------------------------------------------------------

    # DO NOT DELETE OR EDIT NEXT LINE! This is needed for parser.
    # aggregate model marker

    return F
end
