function compute_bcfreq_vardecomp(sr, lr, e_set, m_par; passband=(6,32), ngrid=512)
    indexes_shocks = Vector{Int}(undef, length(e_set.shock_names))
    iter = 1
    SCov = zeros(sr.n_par.nstates_r, sr.n_par.nstates_r)
    for i in e_set.shock_names
        indexes_shocks[iter] = getfield(sr.indexes_r, i)
        SCov[getfield(sr.indexes_r, i), getfield(sr.indexes_r, i)] = (getfield(m_par, Symbol("σ_", i))) .^ 2
        iter += 1
    end
    indexes_states = []
    try
        append!(indexes_states, [sr.indexes_r.distr_m; sr.indexes_r.distr_k; sr.indexes_r.distr_y; sr.indexes_r.COP])
    catch
        append!(indexes_states, [sr.indexes_r.Bhh; sr.indexes_r.Khh])
    end
    for i in sr.state_names
        append!(indexes_states, getfield(sr.indexes_r, Symbol(i)))
    end

    indexes_controls    = setdiff(1:sr.n_par.ntotal_r, indexes_states)
    indexes_endo_states = setdiff(indexes_states, indexes_shocks)
    m = length(indexes_endo_states)
    k = length(indexes_shocks)
    
    # Notation as in Uhlig (2001) (Changes ordering, see below)
    P = lr.LOMstate[indexes_endo_states, indexes_endo_states]
    Q = lr.LOMstate[indexes_endo_states, indexes_shocks]
    R = lr.State2Control[:, indexes_endo_states]
    S = lr.State2Control[:, indexes_shocks]
    N = lr.LOMstate[indexes_shocks, indexes_shocks]

    nvar = sr.n_par.ntotal_r # number of all veriables
    ivar = nvar # number of selected variables (all for now)
 
    freqs = 0:((2*pi)/ngrid):(2*pi*(1-0.5/ngrid)) #[0,2*pi)
    tpos  = exp.(im  .* freqs) #positive frequencies
    tneg  = exp.(-im .* freqs) #negative frequencies

    filter_gain = zeros(1, ngrid)
    lowest_periodicity  = passband[2]
    highest_periodicity = passband[1]
    highest_periodicity = max(2, highest_periodicity) # restrict to upper bound of pi
    filter_gain[(freqs .>=  (2.0 * pi /lowest_periodicity)) .& (freqs .<= (2*pi/highest_periodicity))] .= 1.0
    filter_gain[(freqs .<= (-2.0 * pi /lowest_periodicity + 2 * pi)) .& (freqs .>= (-2*pi/highest_periodicity+2*pi))] .= 1.0


    # Variance decomposition
    VD_alt_order    = zeros(ivar, length(indexes_shocks))

    Sigma_exo = SCov[indexes_shocks, indexes_shocks] + eps()*I# make sure Covariance matrix is positive definite
    
    # Containers for matrices
    mat_bp_col = zeros(ComplexF64, k+1,ngrid, ivar)
    mat1n =zeros(ComplexF64,m+k,k)
    mat1p =zeros(ComplexF64,k,m+k)#[Q' / (I - P' * tpos[1])   I]
    mat2n =zeros(ComplexF64,k,k)
    mat2p =zeros(ComplexF64,k,k)
    mat3n =zeros(ComplexF64,nvar,m+k)
    mat3p =zeros(ComplexF64,m+k,nvar)
    mat31n=zeros(ComplexF64,nvar,k)
    mat13p=zeros(ComplexF64,k,nvar)
    inner =zeros(ComplexF64,k,k)
    f_omega = (1 / (2 * pi)) * inner  # spectral density of state variables; top formula Uhlig [2001], p. 20 with N=0
    g_omega = zeros(ComplexF64,nvar,nvar)
    f_bp = zeros(ComplexF64, ivar, ivar)
    
    for ig = 1:ngrid
        if filter_gain[ig] == 0.0
            mat_bp_col[:,ig, :] = zeros(ComplexF64,k+1, ivar)    # store as matrix row for ifft()
        else
            mat1n .= [(I - P * tneg[ig]) \ Q;  I]
            mat1p .= [Q' / (I - P' * tpos[ig])   I]
            mat2n .= (I - N  * tneg[ig]) 
            mat2p .= (I - N' * tpos[ig]) 
            mat3n .= [I zeros(m, k); R * tneg[ig] S; zeros(k,m) I]
            mat3p .= [I R' * tpos[ig] zeros(m,k); zeros(k,m)  S' I]
            mat31n.= mat3n * mat1n
            mat13p.= mat1p * mat3p
            inner .=((mat2n \ Sigma_exo) / mat2p)
            f_omega .= (1 / (2 * pi)) * inner  # spectral density of state variables; top formula Uhlig [2001], p. 20 with N=0
            g_omega .=  mat31n * f_omega * mat13p # spectral density of selected variables; middle formula Uhlig [2001], p. 20; only middle block, i.e. y_t"
            f_bp  .= filter_gain[ig]^2 * g_omega  # spectral density of selected filtered series; top formula Uhlig [2001], p. 21
            mat_bp_col[1,ig, :] = diag(f_bp)[:]    # store as matrix row for ifft()
            for i=1:k
                Sigma_exo_i      = zeros(size(Sigma_exo))
                Sigma_exo_i[i,i] = Sigma_exo[i,i]
                f_omega         .= (1 / (2 * pi)) * ((mat2n \ Sigma_exo_i) / mat2p)  # spectral density of state variables; top formula Uhlig [2001], p. 20 with N=0
                g_omega         .= mat31n * f_omega * mat13p # spectral density of selected variables; middle formula Uhlig [2001], p. 20; only middle block, i.e. y_t"
                f_bp            .= filter_gain[ig]^2 * g_omega 
                mat_bp_col[i+1,ig, :] = diag(f_bp)[:]
            end
        end    
    end

    imat_bp_col = real.(ifft(mat_bp_col,2)) * (2 * pi) # first dimension shocks (i=1: all), second dimension ngrid, third dimension variables
    UnconditionalVar = imat_bp_col[1, 1, :][:] # Variances
    for i=1:k # Variance contribution of shock i
        VD_alt_order[:,i]= abs.(imat_bp_col[i+1,1, :]) ./ UnconditionalVar 
    end

    # reorder in original format
    reordering        = [indexes_endo_states; indexes_controls; indexes_shocks]
    VD                = similar(VD_alt_order)
    VD[reordering, :] = VD_alt_order
    UnconditionalVar[reordering] = UnconditionalVar

    return VD, UnconditionalVar
end