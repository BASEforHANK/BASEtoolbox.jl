@doc raw"""
    LinearSolution_reduced_system(sr, m_par, A, B,;allow_approx_sol)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys`](@ref), while only differentiating with respect to the
aggregate part of the model, [`Fsys_agg()`](@ref).

The partials of the Jacobian belonging to the heterogeneous agent part of the model
are taken from the full-model derivatives provided as arguments, `A` and `B` (computed
by [`LinearSolution()`](@ref)).

# Arguments
- `sr`: steady-state structure (variable values, indexes, numerical parameters, ...)
- `A`,`B`: derivative of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
- `m_par`: model parameters
- `allow_approx_sol`: if `true`, the function will attempt to solve the linearized model
    even if the system is indeterminate (shifting the critical eigenvalues)

# Returns
as in [`LinearSolution()`](@ref)
"""
function LinearSolution_reduced_system(
    sr,
    m_par,
    A::Array,
    B::Array;
    allow_approx_sol = false,
)

    ############################################################################
    # Calculate dericatives of non-lineear difference equation
    ############################################################################
    length_X0 = length(sr.XSSaggr)
    BA = ForwardDiff.jacobian(
        x -> Fsys_agg(
            x[1:length_X0],
            x[(length_X0 + 1):end],
            sr.XSSaggr,
            sr.distrSS,
            m_par,
            sr.n_par,
            sr.indexes_aggr,
        ),
        zeros(2 * length_X0),
    )

    mask::Vector{Bool} = [!(aggr_names[k] in distr_names) for k in eachindex(aggr_names)]
    indices::Vector{Int} = [getfield(sr.indexes, Symbol(name)) for name in aggr_names[mask]]
    index_XAggP::Vector{Int} = eachindex(aggr_names)[mask] .+ length_X0
    index_XAgg::Vector{Int} = eachindex(aggr_names)[mask]
    # Select only the relevant parts of A,B, and
    grid = [CartesianIndex(i, j) for i in indices for j in indices]
    grid_BA_A = [CartesianIndex(i, j) for i in index_XAgg for j in index_XAggP]
    grid_BA_B = [CartesianIndex(i, j) for i in index_XAgg for j in index_XAgg]
    A[grid] .= BA[grid_BA_A] # copy only the non-distribution part of the Jacobian
    B[grid] .= BA[grid_BA_B] # copy only the non-distribution part of the Jacobian
    # for k in eachindex(aggr_names)
    #     if !(aggr_names[k] in distr_names)
    #         j = indices_k[k]
    #         for h in eachindex(aggr_names)
    #             if !(aggr_names[h] in distr_names)
    #                 i = indices_h[h]
    #                 A[j, i] = BA[k, h+length_X0]
    #                 B[j, i] = BA[k, h]
    #             end
    #         end
    #     end
    # end

    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_LinearSolution, nk = SolveDiffEq(A, B, sr.n_par, allow_approx_sol)

    return gx, hx, alarm_LinearSolution, nk, A, B
end
