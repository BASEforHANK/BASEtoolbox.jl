# function updateV(EVk::Array,
#                  c_a_star::Array,
#                  c_n_star::Array,
#                  m_n_star::Array,
#                  r,
#                  q,
#                  m_par::ModelParameters,
#                  n_par::NumericalParameters
#                 )
#     β::Float64 = m_par.β
#     n = size(c_n_star)
#     #----------------------------------------------------------------------------
#     ## Update Marginal Value Bonds
#     #----------------------------------------------------------------------------
#     Vm = mutil(c_n_star,m_par)                          # marginal utility at consumption policy no adjustment
#     mutil_c_a = mutil(c_a_star,m_par)                          # marginal utility at consumption policy adjustment
#     Vm   .= m_par.λ .* mutil_c_a .+ (1.0 - m_par.λ ) .* Vm # Expected marginal utility at consumption policy (w &w/o adjustment)
    
#     #----------------------------------------------------------------------------
#     ## Update marginal Value of Capital
#     ## i.e. linear interpolate expected future marginal value of capital using savings policy
#     ## Then form expectations.
#     #----------------------------------------------------------------------------
    
#     Vk = similar(EVk)                          # Initialize Vk-container

#     @inbounds @views begin
#         for j::Int = 1:n[3]
#             for k::Int = 1:n[2]
#                  mylinearinterpolate!(Vk[:,k,j],n_par.grid_m, EVk[:,k,j], m_n_star[:,k,j]) # evaluate marginal value at policy
#             end
#         end
#     end
#     Vk        .= r .* Vm .+ m_par.λ .* q .* mutil_c_a .+ (1.0 .- m_par.λ) .* β .* Vk # Expected marginal utility at consumption policy (w &w/o adjustment)
#     return Vk, Vm
# end
function updateV(EVk::Array,
                 c_a_star::Array,
                 c_n_star::Array,
                 m_n_star::Array,
                 r,
                 q,
                 m_par::ModelParameters,
                 n_par::NumericalParameters
                )
    Vk=similar(EVk)
    Vm=similar(EVk)
    mutil_c_a=similar(EVk)
    updateV!(Vk, Vm, mutil_c_a, EVk,
            c_a_star, c_n_star, m_n_star,
            r,  q, m_par, n_par )
    return Vk, Vm
end

function updateV!(Vk::Array,Vm::Array,mutil_c_a::Array,EVk::Array,
    c_a_star::Array, c_n_star::Array, m_n_star::Array,
    r,  q,
    m_par::ModelParameters,
    n_par::NumericalParameters
   )

    β::Float64 = m_par.β
    n = size(c_n_star)
    #----------------------------------------------------------------------------
    ## Update Marginal Value Bonds
    #----------------------------------------------------------------------------
    mutil!(Vm,c_n_star,m_par)                          # marginal utility at consumption policy no adjustment
    mutil!(mutil_c_a,c_a_star,m_par)                          # marginal utility at consumption policy adjustment
    Vm   .= m_par.λ .* mutil_c_a .+ (1.0 - m_par.λ ) .* Vm # Expected marginal utility at consumption policy (w &w/o adjustment)

    #----------------------------------------------------------------------------
    ## Update marginal Value of Capital
    ## i.e. linear interpolate expected future marginal value of capital using savings policy
    ## Then form expectations.
    #----------------------------------------------------------------------------

    #Vk = similar(EVk)                          # Initialize Vk-container

    @inbounds @views begin
        for j::Int = 1:n[3]
            for k::Int = 1:n[2]
                mylinearinterpolate!(Vk[:,k,j],n_par.grid_m, EVk[:,k,j], m_n_star[:,k,j]) # evaluate marginal value at policy
            end
        end
    end
    Vk        .= r .* Vm .+ m_par.λ .* q .* mutil_c_a .+ (1.0 .- m_par.λ) .* β .* Vk # Expected marginal utility at consumption policy (w &w/o adjustment)

end