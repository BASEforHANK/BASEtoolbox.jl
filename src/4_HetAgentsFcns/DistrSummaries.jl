
@doc raw"""
    distrSummaries(distr,c_a_star,c_n_star,n_par,inc,incgross,m_par)

Compute distributional summary statistics, e.g. Gini indexes, top-10%
income and wealth shares, and 10%, 50%, and 90%-consumption quantiles.

# Arguments
- `distr`: joint distribution over bonds, capital and income ``(m \times k \times y)``
- `c_a_star`,`c_n_star`: optimal consumption policies with [`a`] or without [`n`]
    capital adjustment
- `n_par::NumericalParameters`, `m_par::ModelParameters`
- `inc`: vector of (on grid-)incomes, consisting of labor income (scaled by ``\frac{\gamma-\tau^P}{1+\gamma}``, plus labor union-profits),
    rental income, liquid asset income, capital liquidation income,
    labor income (scaled by ``\frac{1-\tau^P}{1+\gamma}``, without labor union-profits),
    and labor income (without scaling or labor union-profits)
- `incgross`: vector of (on grid-) *pre-tax* incomes, consisting of
    labor income (without scaling, plus labor union-profits), rental income,
    liquid asset income, capital liquidation income,
    labor income (without scaling or labor union-profits)
"""
function distrSummaries(distr::AbstractArray,Q::Real,c_a_star::AbstractArray,
                        c_n_star::AbstractArray, n_par::NumericalParameters,
                        inc::AbstractArray,incgross::AbstractArray, m_par::ModelParameters)
    ## Distributional summaries
    distr_m = sum(distr,dims=(2,3))[:]
    distr_k = sum(distr,dims=(1,3))[:]
    distr_y = sum(distr,dims=(1,2))[:]

    total_wealth = Array{eltype(distr)}(undef,n_par.nk.*n_par.nm)
    for k = 1:n_par.nk
        for m = 1:n_par.nm
            total_wealth[m+(k-1)*n_par.nm] = n_par.grid_m[m] .+ Q .* n_par.grid_k[k];
        end
    end
    # Wealth shares and gini
    IX              = sortperm(total_wealth)
    total_wealth    = total_wealth[IX]
    total_wealth_pdf= sum(distr, dims=3)
    total_wealth_pdf= total_wealth_pdf[IX];
    total_wealth_cdf= cumsum(total_wealth_pdf);
    total_wealth_w  = total_wealth.*total_wealth_pdf # weighted
    wealthshares    = cumsum(total_wealth_w)./sum(total_wealth_w);

    TOP10Wshare     = 1.0 - mylinearinterpolate(total_wealth_cdf,wealthshares,[0.9])[1]
    giniwealth      = gini(total_wealth,total_wealth_pdf)
    # Consumption distribution
    c               = Array{eltype(c_a_star)}(undef,(n_par.nm, n_par.nk,n_par.ny,2))
    distr_c         = similar(c)
    aux_x           = inc[5] # adjustment for labor in GHH preferences
    aux_x[:,:,end]  .= zeros(n_par.nm,n_par.nk); # entrepreneurs do not work
    c[:,:,:,1]      .= c_a_star .+ aux_x; # add adjustment to consumption
    c[:,:,:,2]      .= c_n_star .+ aux_x; # add adjustment to consumption
    distr_c[:,:,:,1].= m_par.λ .* distr
    distr_c[:,:,:,2].= (1-m_par.λ) .* distr

    # Gini of goods consumption
    IX              = sortperm(c[:]);
    c[:]            .= c[IX];
    distr_c[:]      .= distr_c[IX];
    #S               = [0.0; cumsum(c_pdf.*c)];
    #giniconsumption = 1.0 .- (sum(c_pdf.*(S[1:end-1] .+ S[2:end]))/S[end]);
    giniconsumption = gini(c,distr_c)

    # Top 10 net income share
    capital_inc     = inc[2] .+ inc[3] .- n_par.mesh_m
    Yidio           = inc[6] .+ capital_inc 
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    Y_pdf           = distr[IX]
    Y_cdf           = cumsum(Y_pdf)
    Y_w             = Yidio.*Y_pdf
    net_incomeshares= cumsum(Y_w)./sum(Y_w);
    TOP10Inetshare  = 1.0 .- mylinearinterpolate(Y_cdf,net_incomeshares,[0.9])[1]

    # Top 10 gross income share
    Yidio           = incgross[1] .+ capital_inc 
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    Y_pdf           = distr[IX]
    Y_cdf           = cumsum(Y_pdf)
    Y_w             = Yidio.*Y_pdf
    incomeshares    = cumsum(Y_w)./sum(Y_w);
    TOP10Ishare     = 1.0 .- mylinearinterpolate(Y_cdf,incomeshares,[0.9])[1]

    # Standard deviation of log labor earnings
    Yidio           = log.(incgross[1][:,:,1:end-1])
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    distr_aux       = distr[:,:,1:end-1]
    distr_aux       = distr_aux./sum(distr_aux[:])
    Y_pdf           = distr_aux[IX]
   
    sdlogy          = sqrt(dot(Y_pdf,Yidio.^2) .- dot(Y_pdf, Yidio).^2);



    return     distr_m, distr_k, distr_y,  TOP10Wshare, TOP10Ishare,TOP10Inetshare, 
               giniwealth, giniconsumption, sdlogy
end
function gini(x,pdf)
    s=0.0
    gini=0.0
    for i =eachindex(x)
        gini -= pdf[i]*s
        s    += x[i]*pdf[i]
        gini -= pdf[i]*s
    end
    gini /=s
    gini +=1.0
    return gini
end