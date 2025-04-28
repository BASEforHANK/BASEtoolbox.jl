
@doc raw"""
    model_reduction()

Produce Model Reduction based on Variance Covariance Matrix of States and Controls.

# Returns/ Updates
`SteadyResults`, containing returns of [`find_steadystate()`](@ref)
"""
function model_reduction(sr, lr, m_par)
    @printf "\n"
    @printf "Model reduction (state-space representation)...\n"

    n_par = sr.n_par
    # Reduce further based on importance in dynamics at initial guess
    if n_par.further_compress
        @printf "Reduction Step\n"
        indexes_r, n_par = compute_reduction(sr, lr, m_par, e_set.shock_names)

        @printf "Number of reduced model factors for DCTs for Wb & Wk: %d\n" (
            length(indexes_r.Wb) + length(indexes_r.Wk)
        )

        @printf "Number of reduced model factors for copula DCTs: %d\n" length(
            indexes_r.COP,
        )
    else
        @printf "Further model reduction switched off --> reverting to full model\n"
        @set! n_par.PRightAll = Diagonal(ones(n_par.ntotal))#float(I[1:n_par.ntotal, 1:n_par.ntotal])
        @set! n_par.PRightStates = Diagonal(ones(n_par.nstates))# float(I[1:n_par.nstates, 1:n_par.nstates])
        indexes_r = sr.indexes
        @set! n_par.nstates_r = n_par.nstates
        @set! n_par.ncontrols_r = n_par.ncontrols
        @set! n_par.ntotal_r = n_par.ntotal
    end

    @printf "Model reduction (state-space representation)... Done.\n"

    return SteadyResults(
        sr.XSS,
        sr.XSSaggr,
        sr.indexes,
        indexes_r,
        sr.indexes_aggr,
        sr.compressionIndexes,
        n_par,
        m_par,
        sr.CDFSS,
        sr.CDF_bSS,
        sr.CDF_kSS,
        sr.CDF_hSS,
        sr.distrSS,
        state_names,
        control_names,
    )
end
