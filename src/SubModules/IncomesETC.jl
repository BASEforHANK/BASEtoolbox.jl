module IncomesETC

using ..Parsing

using LinearAlgebra, Roots, ForwardDiff, Printf

export incomes!,
    incomes,
    av_labor_tax_rate,
    mutil!,
    mutil,
    invmutil!,
    invmutil,
    employment,
    interest,
    wage,
    output,
    profits_U,
    CompMarketsCapital,
    labor_supply,
    compute_args_hh_prob_ss,
    scale_Hprog,
    borrowing_rate_ss

# Scaling of labor disutility including tax progressivity, see eq. (Gross income)
scale_Hprog(τprog, m_par) = (m_par.γ .+ (m_par.Tprog .- 1.0)) ./ (m_par.γ .+ τprog)

# Scaling for composite good, see eq. (BC with x)
scale_GHH(τprog, m_par) = (m_par.γ .+ τprog) ./ (1.0 .+ m_par.γ)

# Tax functions, see eq. (Tax func) and text below
labor_tax_f(x, τlev, τprog, tax_base, scaling::Bool) =
    (1.0 .- τlev) .* x .^ (1.0 .- τprog) .* (tax_base^τprog) .^ scaling
union_tax_f(x, τbar) = (1.0 .- τbar) .* x

# Labor supply function, see eq. (Total hours)
function labor_supply(wH, Hprog, τlev, τprog, τc, m_par, tax_base, scaling::Bool)
    if scaling == true && tax_base == 0.0
        @warn "Scaling switched on but tax base is zero. Are you sure?"
    end

    return (
        ((1.0 .- τprog) .* (1.0 .- τlev) ./ (1.0 .+ τc) .* (tax_base^τprog) .^ scaling) .^
        (1.0 ./ (m_par.γ .+ τprog)) .* (wH) .^ ((1.0 .- τprog) ./ (m_par.γ .+ τprog)) .*
        Hprog
    )
end

include("IncomesETC/incomes.jl")
include("IncomesETC/av_labor_tax_rate.jl")
include("IncomesETC/other.jl")

# Documentation mode: if paths to model are not defined, the code will use the baseline example.
if !isdefined(Main, :paths)
    include("../../examples/baseline/Model/input_functions.jl")
    include("../../examples/baseline/Model/input_compute_args_hh_prob_ss.jl")
else
    include(Main.paths["src_example"] .* "/Model/input_functions.jl")
    include(Main.paths["src_example"] .* "/Model/input_compute_args_hh_prob_ss.jl")
end

end
