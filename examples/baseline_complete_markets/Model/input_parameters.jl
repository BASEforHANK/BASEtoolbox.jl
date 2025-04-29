@metadata prior nothing
@metadata label ""
@metadata latex_label L""
@metadata long_name ""

"""
    ModelParameters()

A structure to collect all model parameters, including calibrated values and prior
distributions for estimation.

# Overview

  - This struct is designed for macroeconomic models and includes parameters related to
    household preferences, income processes, technological factors, monetary policy, fiscal
    policy, and exogenous shocks.
  - The parameters are annotated with metadata such as names (both ASCII and LaTeX), prior
    distributions, and a boolean flag indicating whether they are estimated.
  - Uses the `Parameters`, `FieldMetadata`, and `Flatten` packages to facilitate parameter
    management.

# Fields

Each field follows the structure:

```julia
parameter::T = value | "ascii_name" | L"latex_name" | prior_distribution | estimated
```

  - `parameter`: Internal model parameter name.
  - `value`: Default numerical value.
  - `ascii_name`: Human-readable name used in output or logging.
  - `latex_name`: Corresponding LaTeX notation for use in reports and documentation.
  - `prior_distribution`: Prior distribution for Bayesian estimation (if applicable).
  - `estimated`: Boolean indicating whether the parameter is estimated.
"""
@label @long_name @latex_label @prior @flattenable @with_kw struct ModelParameters{T}

    # Household preference parameters
    ξ::T = 4.0 | "xi" | "risk aversion" | L"\xi" | _ | false
    γ::T = 2.0 | "gamma" | "inverse Frisch elasticity" | L"\gamma" | _ | false
    β::T = 0.98255 | "beta" | "discount factor" | L"\beta" | _ | false
    λ::T =
        0.065 | "lambda" | "illiquid asset adjustment probability" | L"\lambda" | _ | false

    # Individual income process
    ρ_h::T = 0.98 | "rho" | "autocorrelation of income shocks" | L"\rho" | _ | false
    σ_h::T = 0.12 | "sigma" | "standard deviation of income shocks" | L"\sigma" | _ | false
    ι::T = 1 / 16 | "iota" | "probability to return worker" | L"\iota" | _ | false
    ζ::T = 1 / 4500 | "zeta" | "probability to become entrepreneur" | L"\zeta" | _ | false

    # Technological parameters
    α::T = 0.318 | "alpha" | "capital share" | L"\alpha" | _ | false
    δ_0::T = 0.0215 | "delta" | "depreciation rate" | L"\delta" | _ | false
    δ_s::T = 0.70557 | "delta_s" | "slope of depreciation rate" | L"\delta_s" | _ | false
    ϕ::T = 1.94092 | "phi" | "capital adjustment costs" | L"\phi" | _ | false
    μ::T = 1.1 | "mu" | "price markup" | L"\mu" | _ | false
    κ::T = 0.14561 | "kappa" | "price adjustment costs" | L"\kappa" | _ | false
    μw::T = 1.1 | "mu_w" | "wage markup" | L"\mu_w" | _ | false
    κw::T = 0.23931 | "kappa_w" | "wage adjustment costs" | L"\kappa_w" | _ | false

    # Further steady-state parameters
    Tlev::T = 1.1775 | "tau_lev" | "income tax rate level" | L"\tau^l" | _ | false
    Tprog::T = 1.1022 | "tau_pro" | "income tax rate progressivity" | L"\tau^p" | _ | false
    Tc::T = 1.0 | "Tc" | "VAT rate (gross)" | L"T_c" | _ | false
    RRB::T = 1.0 | "RB" | "real rate on bonds (gross)" | L"RRB" | _ | false
    Rbar::T = 0.02178 | "Rbar" | "borrowing wedge" | L"\bar R" | _ | false
    q::T = 1.0 | "q" | "price of capital" | L"q" | _ | false
    Z::T = 1.0 | "Z" | "TFP" | L"Z" | _ | false
    σ::T = 1.0 | "sigma" | "income risk" | L"\sigma" | _ | false
    ψ::T = 0.25 | "psi" | "share of capital in liquid assets" | L"\psi" | _ | false

    # Tradable shares
    ωΠ::T = 0.2 | "omegaPi" | "fraction tradable profits" | L"\omega^{\Pi}" | _ | false
    ιΠ::T = 0.016 | "iotaPi" | "fraction depreciating shares" | L"\iota^{\Pi}" | _ | false

    # monetary policy
    ρ_R::T = 0.80306 | "rho_R" | "persistence of Taylor rule" | L"\rho_R" | _ | false
    θ_π::T =
        2.07808 | "theta_pi" | "reaction of Taylor rule to inflation" | L"\theta_\pi" | _ |
        false
    θ_Y::T =
        0.21873 | "theta_Y" | "reaction of Taylor rule to output" | L"\theta_y" |
        Normal(0.125, 0.05) | false

    # fiscal policy
    scale_prog::Bool =
        false | "scale_prog" | "scaling of tax rate with tax base" | "scale_prog" | _ |
        false
    γ_B::T = 0.102013 | "gamma_B" | "reaction of deficit to debt" | L"\gamma_B" | _ | false
    γ_π::T =
        -2.1737 | "gamma_pi" | "reaction of deficit to inflation" | L"\gamma_{\pi}" | _ |
        false
    γ_Y::T =
        -0.43631 | "gamma_Y" | "reaction of deficit to output" | L"\gamma_Y" | _ | false

    ρ_τ::T = 0.49265 | "rho_tau" | "persistence of tax rule" | L"\rho_\tau" | _ | false
    γ_Bτ::T =
        3.29306 | "gamma_Btau" | "reaction of tax level to debt" | L"\gamma_B^\tau" | _ |
        false
    γ_Yτ::T =
        -0.92073 | "gamma_Ytau" | "reaction of tax level to output" | L"\gamma_Y_\tau" | _ |
        false

    ρ_P::T = 0.91942 | "rho_P" | "persistence of tax progressivity" | L"\rho_P" | _ | false
    γ_BP::T =
        0.0 | "gamma_BP" | "reaction of tax progressivity to debt" | L"\gamma_B^P" | _ |
        false
    γ_YP::T =
        0.0 | "gamma_YP" | "reaction of tax progressivity to output" | L"\gamma_Y^P" | _ |
        false
    γ_WP::T =
        0.0 | "gamma_WP" | "reaction of tax progressivity to output" | L"\gamma_W^P" | _ |
        false

    # exogeneous aggregate "shocks"

    ρ_Sshock::T =
        1e-8 | "rho_Sshock" | "autocorrelation of uncertainty shock" | L"\rho_{Sshock}" |
        _ | false

    ρ_Tprogshock::T =
        1e-8 | "rho_Pshock" | "autocorrelation of tax progressivity shock" |
        L"\rho_{Pshock}" | _ | false
    σ_Tprogshock::T =
        0.06865 | "sigma_Pshock" | "standard deviation of tax progressivity shock" |
        L"\sigma_P" | _ | false

    ρ_Gshock::T =
        0.96822 | "rho_Gshock" | "autocorrelation of deficit shock" | L"\rho_D" | _ | false
    σ_Gshock::T =
        0.00376 | "sigma_G" | "standard deviation of deficit shock" | L"\sigma_D" | _ |
        false

    ρ_Rshock::T =
        1e-8 | "rho_Rshock" | "autocorrelation of Taylor rule shock" | L"\rho_{Rshock}" |
        _ | false
    σ_Rshock::T =
        0.00231 | "sigma_Rshock" | "standard deviation of Taylor rule shock" | L"\sigma_R" |
        _ | false

    ρ_A::T =
        0.97241 | "rho_A" | "autocorrelation of bond-spread shock" | L"\rho_A" | _ | false
    σ_A::T =
        0.00158 | "sigma_A" | "standard deviation of bond-spread shock" | L"\sigma_A" | _ |
        false

    ρ_Z::T = 0.99782 | "rho_Z" | "autocorrelation of TFP shock" | L"\rho_Z" | _ | false
    σ_Z::T =
        0.00601 | "sigma_Z" | "standard deviation of TFP shock" | L"\sigma_Z" | _ | false

    ρ_ZI::T =
        0.76371 | "rho_Psi" | "autocorrelation of investment efficiency shock" |
        L"\rho_\Psi" | _ | false
    σ_ZI::T =
        0.07211 | "sigma_Psi" | "standard deviation of investment efficiency shock" |
        L"\sigma_\Psi" | _ | false

    ρ_μ::T =
        0.90374 | "rho_mu" | "autocorrelation of price markup shock" | L"\rho_\mu" | _ |
        false
    σ_μ::T =
        0.01351 | "sigma_mu" | "standard deviation of price markup shock" | L"\sigma_\mu" |
        _ | false

    ρ_μw::T =
        0.90579 | "rho_muw" | "autocorrelation of wage markup shock" | L"\rho_{\mu w}" | _ |
        false
    σ_μw::T =
        0.03506 | "sigma_muw" | "standard deviation of wage markup shock" |
        L"\sigma_{\mu w}" | _ | false

    ρ_s::T =
        0.54472 | "rho_sigma" | "autocorrelation of income risk shock" | L"\rho_s" | _ |
        false
    σ_Sshock::T =
        0.69186 | "sigma_Sshock" | "standard deviation of income risk shock" | L"\sigma_s" |
        _ | false
    Σ_n::T =
        28.87977 | "Sigma_n" | "reaction of income risk to aggregates" | L"\Sigma_N" | _ |
        false
end

"""
    NumericalParameters()

Collect parameters for the numerical solution of the model in a `struct`.
"""
@with_kw struct NumericalParameters
    m_par::ModelParameters = ModelParameters()

    # If there are crashes in linear interpolations, then this often is because the grid is
    # not ideal or the interest rate is below 0. This can happen when the algorithm searches
    # for an equilibrium but does not need to be a serious mistake. To avoid the crashes,
    # set the following parameter to true, but consider adjusting the maximal gris
    # parameter.
    warn_egm::Bool = true

    # model we are solving
    model = CompleteMarkets()

    # regular grid
    nh::Int = 2
    nk::Int = 1
    nb::Int = 1

    # copula grid
    nh_copula::Int = 2 # rule of thumb: divide nh, without entrepreneur, by two
    nk_copula::Int = 1 # rule of thumb: divide nk by twelve
    nb_copula::Int = 1 # rule of thumb: divide nb by twelve

    # coarse grid in find_steadystate
    nh_coarse::Int = 2
    nk_coarse::Int = 1
    nb_coarse::Int = 1

    # capital bounds for coarse grid in find_steadystate
    Kmin_coarse::Float64 = 20.0
    Kmax_coarse::Float64 = 51.0

    # other options for find_steadystate
    rKmin_coarse::Float64 = 0.0
    search_range::Float64 = 0.1

    # minimum and maximum values for grids
    kmin::Float64 = typeof(model) == TwoAsset ? 0.0 : 0.0
    kmax::Float64 = typeof(model) == TwoAsset ? 4500.0 : 0.0
    bmin::Float64 = typeof(model) == CompleteMarkets ? 0.0 : -6.0
    bmax::Float64 = typeof(model) == CompleteMarkets ? 0.0 : 3500.0

    # precision of solution
    ϵ::Float64 = 1e-13

    sol_algo::Symbol = :schur # options: :schur (Klein's method), :lit (linear time iteration), :litx (linear time iteration with Howard improvement)
    verbose::Bool = true   # verbose model
    reduc_value::Float64 = 1e-3   # Lost fraction of "energy" in the DCT compression for value functions
    reduc_marginal_value::Float64 = 1e-3   # Lost fraction of "energy" in the DCT compression for value functions

    further_compress::Bool = true   # run model-reduction step based on MA(∞) representation
    further_compress_critC = eps()  # critical value for eigenvalues for Value functions
    further_compress_critS = ϵ      # critical value for eigenvalues for copula

    # Parameters that will be overwritten in the code
    aggr_names::Array{String,1} = ["Something"] # Placeholder for names of aggregates
    distr_names::Array{String,1} = ["Something"] # Placeholder for names of distributions

    naggrstates::Int = 16 # (placeholder for the) number of aggregate states
    naggrcontrols::Int = 16 # (placeholder for the) number of aggregate controls
    nstates::Int = nh + nk + nb + naggrstates - 3 # (placeholder for the) number of states + controls in total
    ncontrols::Int = 16 # (placeholder for the) number of controls in total
    ntotal::Int = nstates + ncontrols     # (placeholder for the) number of states+ controls in total
    n_agg_eqn::Int = nstates + ncontrols     # (placeholder for the) number of aggregate equations
    naggr::Int = length(aggr_names)     # (placeholder for the) number of aggregate states + controls
    ntotal_r::Int = nstates + ncontrols# (placeholder for the) number of states + controls in total after reduction
    nstates_r::Int = nstates# (placeholder for the) number of states in total after reduction
    ncontrols_r::Int = ncontrols# (placeholder for the) number of controls in total after reduction

    PRightStates::AbstractMatrix = Diagonal(ones(nstates)) # (placeholder for the) Matrix used for second stage reduction (states only)
    PRightAll::AbstractMatrix = Diagonal(ones(ntotal))  # (placeholder for the) Matrix used for second stage reduction

    # Transition matrix and grid for income in steady state (worker - entrepreneur)
    grid_h::Array{Float64,1} = [
        1.0
        (m_par.ζ .+ m_par.ι) / m_par.ζ
    ]
    Π::Matrix{Float64} = [
        (1.0 .- m_par.ζ) m_par.ζ
        m_par.ι 1.0 .- m_par.ι
    ]
    # bounds of income bins (except entrepreneur)
    bounds_h::Array{Float64,1} = Tauchen(m_par.ρ_h, nh - 1)[3]

    Htilde::Float64 = ((Π ^ 1000)[1, 1:(end - 1)]' * grid_h[1:(end - 1)]) # stationary equilibrium average human capital
    frac_workers::Float64 = (1.0 / (1.0 - (Π ^ 1000)[1, end]))     # stationary equilibrium fraction workers

    # initial gues for stationary distribution (needed if iterative procedure is used)
    dist_guess::Array{Float64,3} = ones(nb, nk, nh) / (nb * nk * nh)

    # grid illiquid assets:
    grid_k::Array{Float64,1} =
        exp.(range(log(kmin + 1.0); stop = log(kmax + 1.0), length = nk)) .- 1.0

    # grid liquid assets:
    grid_b::Array{Float64,1} =
        exp.(range(0; stop = log(bmax - bmin + 1.0), length = nb)) .+ bmin .- 1.0

    # meshes for income, liquid and illiquid assets
    mesh_h::Array{Float64,3} = repeat(reshape(grid_h, (1, 1, nh)); outer = [nb, nk, 1])
    mesh_b::Array{Float64,3} = repeat(reshape(grid_b, (nb, 1, 1)); outer = [1, nk, nh])
    mesh_k::Array{Float64,3} = repeat(reshape(grid_k, (1, nk, 1)); outer = [nb, 1, nh])

    # grid for copula marginal distributions
    copula_marginal_b::Array{Float64,1} =
        nb == 1 ? [1.0] : collect(range(0.0; stop = 1.0, length = nb_copula))
    copula_marginal_k::Array{Float64,1} =
        nk == 1 ? [1.0] : collect(range(0.0; stop = 1.0, length = nk_copula))
    copula_marginal_h::Array{Float64,1} =
        nh == 1 ? [1.0] : collect(range(0.0; stop = 1.0, length = nh_copula))

    # Storage for linearization results
    LOMstate_save::Array{Float64,2} = zeros(nstates, nstates)
    State2Control_save::Array{Float64,2} = zeros(ncontrols, nstates)
end

"""
    EstimationSettings()

Collect settings for the estimation of the model parameters in a `struct`.

Use package `Parameters` to provide initial values. Input and output file names are stored
in the fields `mode_start_file`, `data_file`, `save_mode_file` and `save_posterior_file`.
"""
@with_kw struct EstimationSettings
    shock_names::Array{Symbol,1} = shock_names # set in Model/input_aggregate_names.jl
    observed_vars_input::Array{Symbol,1} = [
        :Ygrowth,
        :Igrowth,
        :Cgrowth,
        :N,
        :wgrowth,
        :RB,
        :π,
        :TOP10Wshare,
        :TOP10Ishare,
        :τprog,
        :σ,
    ]

    nobservables = length(observed_vars_input)

    data_rename::Dict{Symbol,Symbol} = Dict(
        :pi => :π,
        :sigma2 => :σ,
        :tauprog => :τprog,
        :w90share => :TOP10Wshare,
        :I90share => :TOP10Ishare,
    )

    me_treatment::Symbol = :unbounded
    me_std_cutoff::Float64 = 0.2

    meas_error_input::Array{Symbol,1} = [:TOP10Wshare, :TOP10Ishare]
    meas_error_distr::Array{InverseGamma{Float64},1} =
        [InverseGamma(ig_pars(0.01, 0.01^2)...), InverseGamma(ig_pars(0.01, 0.01^2)...)]

    # Leave empty to start with prior mode
    mode_start_file::String = ""

    data_file::String = ""
    save_mode_file::String = ""
    save_posterior_file::String = ""

    estimate_model::Bool = true

    max_iter_mode::Int = 3
    optimizer::Optim.AbstractOptimizer = NelderMead()
    compute_hessian::Bool = false    # true: computes Hessian at posterior mode; false: sets Hessian to identity matrix
    f_tol::Float64 = 1.0e-4
    x_tol::Float64 = 1.0e-4

    multi_chain_init::Bool = false
    ndraws::Int = 400
    burnin::Int = 100
    mhscale::Float64 = 0.00015
    debug_print::Bool = true
    seed::Int = 778187
end
