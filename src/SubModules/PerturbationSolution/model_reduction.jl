
@doc raw"""
    model_reduction()

Produce Model Reduction based on Variance Covariance Matrix of States and Controls.

# Returns/ Updates
`struct` `SteadyResults`, containing returns of [`find_steadystate()`](@ref)
"""
function model_reduction(sr, lr, m_par)
    n_par = sr.n_par
    # Reduce further based on importance in dynamics at initial guess 
    if n_par.further_compress
        println("Reduction Step")
        indexes_r, n_par = compute_reduction(sr, lr, m_par, e_set.shock_names)

        println("Number of reduced model factors for DCTs for Vm & Vk:")
        println(length(indexes_r.Vm) + length(indexes_r.Vk))

        println("Number of reduced model factors for copula DCTs:")
        println(length(indexes_r.COP))
    else
        println("Further model reduction switched off --> reverting to full model")
        @set! n_par.PRightAll = Diagonal(ones(n_par.ntotal))#float(I[1:n_par.ntotal, 1:n_par.ntotal])
        @set! n_par.PRightStates = Diagonal(ones(n_par.nstates))# float(I[1:n_par.nstates, 1:n_par.nstates])
        indexes_r = sr.indexes
        @set! n_par.nstates_r = n_par.nstates
        @set! n_par.ncontrols_r = n_par.ncontrols
        @set! n_par.ntotal_r = n_par.ntotal
    end

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
        sr.CDF_m,
        sr.CDF_k,
        sr.CDF_y,
        sr.distrSS,
        state_names,
        control_names,
    )
end

