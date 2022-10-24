function incomes(n_par, m_par, mcw, A, q, RB, τprog, τlev, H, Ht, π,r,w,N, profits, unionprofits, av_tax_rate)
# consider ! version that does not allocate inc new
    inc = fill(Array{typeof(RB)}(undef, size(n_par.mesh_m)),6)
    incgross = fill(Array{typeof(RB)}(undef, size(n_par.mesh_m)),5)
    eff_int = Array{typeof(RB)}(undef, size(n_par.mesh_m))
    incomes!(inc,incgross,eff_int, n_par, m_par, mcw, A, q, RB, τprog, τlev, H, Ht, π,r,w,N, profits, unionprofits, av_tax_rate)
return incgross, inc, eff_int
end
function incomes!(inc, incgross, eff_int, n_par, m_par, mcw, A, q, RB, τprog, τlev, H, Ht, π,r,w,N, profits, unionprofits, av_tax_rate)
    tax_prog_scale = (m_par.γ + m_par.τ_prog)/((m_par.γ + τprog))
    labor_inc      = mcw.*w.*N./Ht.*(n_par.mesh_y/H).^tax_prog_scale 
    net_labor_inc  = τlev .* labor_inc.^(1.0-τprog)
    net_u_profits  = unionprofits .* (1.0 .- av_tax_rate).* n_par.HW
    eff_int        .= ((RB .* A) ./ π        .+ (m_par.Rbar .* (n_par.mesh_m.<=0.0)))  # effective rate (need to check timing below and inflation)
    GHHFA = ((m_par.γ + τprog)/(m_par.γ+1)) # transformation (scaling) for composite good

    inc   .= [   GHHFA.* net_labor_inc .+ net_u_profits, # incomes of workers adjusted for disutility of labor
                (r .- 1.0)   .* n_par.mesh_k,                      # rental income
                eff_int      .* n_par.mesh_m,                      # liquid asset Income
                q            .* n_par.mesh_k,                      # capital liquidation Income (q=1 in steady state)
                τlev .* ((1.0 - τprog)/(m_par.γ+1)) .* (mcw.*w.*N.*n_par.mesh_y./ H).^(1.0-τprog),
                net_labor_inc .+ net_u_profits
                ] 
    inc[1][:,:,end].= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes
    inc[5][:,:,end].= 0.0
    inc[6][:,:,end].= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes

    incgross .= [ labor_inc .+ unionprofits .* n_par.HW,    # labor income
                inc[2],                 # rental income
                inc[3],                 # liquid asset Income
                inc[4],                 # capital liquidation Income (q=1 in steady state)
                labor_inc               # labor income w/o union
                ]           
    incgross[1][:,:,end].= (n_par.mesh_y[:,:,end] .* profits)
    incgross[5][:,:,end].= (n_par.mesh_y[:,:,end] .* profits)

#return incgross, inc, eff_int
end