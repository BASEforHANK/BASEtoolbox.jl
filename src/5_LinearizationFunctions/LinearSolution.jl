@doc raw"""
    LinearSolution(sr, m_par, A, B; estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys()`](@ref), using Schmitt-Grohé & Uribe (JEDC 2004) style linearization
(apply the implicit function theorem to obtain linear observation and
state transition equations).

The Jacobian is calculated using the package `ForwardDiff`

# Arguments
- `sr`: steady-state structure (variable values, indexes, numerical parameters, ...)
- `A`,`B`: derivative of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
- `m_par`: model parameters

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_LinearSolution`,`nk`: `alarm_LinearSolution=true` when solving algorithm fails, `nk` number of
    predetermined variables
- `A`,`B`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
"""
function LinearSolution(sr::SteadyResults, m_par::ModelParameters, A::Array, B::Array; estim=false)
    ############################################################################
    # Prepare elements used for uncompression
    ############################################################################
    # Matrices to take care of reduced degree of freedom in marginal distributions
    Γ  = shuffleMatrix(sr.distrSS, sr.n_par)
    # Matrices for discrete cosine transforms
    DC = Array{Array{Float64,2},1}(undef,3)
    DC[1]  = mydctmx(sr.n_par.nm)
    DC[2]  = mydctmx(sr.n_par.nk)
    DC[3]  = mydctmx(sr.n_par.ny)
    IDC    = [DC[1]', DC[2]', DC[3]']

    DCD = Array{Array{Float64,2},1}(undef,3)
    DCD[1]  = mydctmx(sr.n_par.nm_copula)
    DCD[2]  = mydctmx(sr.n_par.nk_copula)
    DCD[3]  = mydctmx(sr.n_par.ny_copula)
    IDCD    = [DCD[1]', DCD[2]', DCD[3]']

    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################
    length_X0 = sr.n_par.ntotal 
    X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,0.0)
    F  = Fsys(X0, X0, sr.XSS, m_par, sr.n_par, sr.indexes, Γ, sr.compressionIndexes, DC, IDC, DCD, IDCD)
   
    FR=realpart.(F)
    println(findall(abs.(FR).>0.001))
    println("Number of States and Controls")
    println(length(F))
    println("Max error on Fsys:")
    println(maximum(abs.(FR[:])))
    println("Max error of COP in Fsys:")
    println(maximum(abs.(FR[sr.indexes.COP])))
    println("Max error of Vm in Fsys:")
    println(maximum(abs.(FR[sr.indexes.Vm])))

    println("Max error of Vk in Fsys:")
    println(maximum(abs.(FR[sr.indexes.Vk])))
    
    ############################################################################
    # Calculate Jacobians of the Difference equation F
    ############################################################################
    # BA  = ForwardDiff.jacobian(x-> Fsys(x[1:length_X0], x[length_X0+1:end],
    #                 sr.XSS, m_par, sr.n_par, sr.indexes, Γ, sr.compressionIndexes, DC, IDC, DCD, IDCD), zeros(2*length_X0))

    # B   = BA[:,1:length_X0]
    # A   = BA[:,length_X0+1:end]
    f(x) = Fsys(x[1:length_X0], x[length_X0+1:end], sr.XSS, m_par, 
                sr.n_par, sr.indexes, Γ, sr.compressionIndexes, DC, IDC, DCD, IDCD)
    
    function manip_some(x,indexes, lengthy)
        y = zeros(eltype(x),lengthy) 
        y[indexes] = x
        return y
    end
    A = zeros(length_X0,length_X0)
    B = zeros(length_X0,length_X0)

    dist_indexes    = [sr.indexes.distr_m; sr.indexes.distr_k; sr.indexes.distr_y; sr.indexes.COP] # changes in marginal distributions at time t+1 affect t+1 copula errors
    V_indexes       = [sr.indexes.Vm; sr.indexes.Vk]
    
    not_dist_indexes= setdiff(1:length_X0,dist_indexes)
    not_V_indexes   = setdiff(1:length_X0,V_indexes)
    # Derivatives with respect to time t+1 distributions are known (unit/shuffle mat)
    A[dist_indexes , dist_indexes] = -I[1:length(dist_indexes), 1:length(dist_indexes)]
    A[sr.indexes.distr_m,sr.indexes.distr_m] = -Γ[1][1:end-1,:]
    A[sr.indexes.distr_k,sr.indexes.distr_k] = -Γ[2][1:end-1,:]
    A[sr.indexes.distr_y,sr.indexes.distr_y] = -Γ[3][1:end-1,:]
    # Derivatives w.r.t. time t value functions are known (unit matrix)
    B[V_indexes    , V_indexes]    = I[1:length(V_indexes)   , 1:length(V_indexes)]
    
    A[:, not_dist_indexes]  = ForwardDiff.jacobian(x-> f(manip_some(x,not_dist_indexes .+ length_X0, 2*length_X0)), zeros(length(not_dist_indexes)))
    B[:, not_V_indexes]     = ForwardDiff.jacobian(x-> f(manip_some(x,not_V_indexes, 2*length_X0)), zeros(length(not_V_indexes)))

    ######################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_LinearSolution, nk = SolveDiffEq(A, B, sr.n_par, estim)

    println("State Space Solution Done")

    return gx, hx, alarm_LinearSolution, nk, A, B
end

