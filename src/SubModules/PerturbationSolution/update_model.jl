
@doc raw"""
    update_model()

Updates the linearized model (around the steady state, after parameter changes in the
aggregate model) and solves, using [`LinearSolution_reduced_system()`](@ref). WARNING: The
function is not threadsafe in the sense that calling it will alter the input(!) lr.A/B
across threads, if lr is not local to the thread.

# Returns
`LinearResults`, containing
- `A::Array{Float64,2}`,`B::Array{Float64,2}`: first derivatives of [`Fsys()`](@ref) with
    respect to arguments `X` [`B`] and `XPrime` [`A`]
- `State2Control::Array{Float64,2}`: observation equation
- `LOMstate::Array{Float64,2}`: state transition equation
"""
function update_model(sr::SteadyResults, lr::LinearResults, m_par::ModelParameters)
    if sr.n_par.verbose
        @printf "Updating linearization\n"
    end
    State2Control, LOMstate, SolutionError, nk, A, B =
        LinearSolution_reduced_system(sr, m_par, lr.A, lr.B; allow_approx_sol = false)

    return LinearResults(State2Control, LOMstate, A, B, SolutionError, nk)
end
