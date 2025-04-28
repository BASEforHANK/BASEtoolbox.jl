"""
Mainboard for the baseline example of the BASEforHANK package.
"""

using PrettyTables, Printf, BenchmarkTools;

## ------------------------------------------------------------------------------------------
## Header: set up paths, pre-process user inputs, load module
## ------------------------------------------------------------------------------------------

root_dir = replace(Base.current_project(), "Project.toml" => "");
cd(root_dir);

# set up paths for the project
paths = Dict(
    "root" => root_dir,
    "src" => joinpath(root_dir, "src"),
    "bld" => joinpath(root_dir, "bld"),
    "src_example" => @__DIR__,
    "bld_example" => replace(@__DIR__, "examples" => "bld") * "_estim",
);

# create bld directory for the current example
mkpath(paths["bld_example"]);

# pre-process user inputs for model setup
include(paths["src"] * "/Preprocessor/PreprocessInputs.jl");
include(paths["src"] * "/BASEforHANK.jl");
using .BASEforHANK;

# set BLAS threads to the number of Julia threads, prevents grabbing all
BASEforHANK.LinearAlgebra.BLAS.set_num_threads(Threads.nthreads());

## ------------------------------------------------------------------------------------------
## Initialize: set up model parameters, priors, and estimation settings
## ------------------------------------------------------------------------------------------

# model parameters and priors
m_par = ModelParameters();
priors = collect(metaflatten(m_par, prior));
par_prior = mode.(priors);
m_par = BASEforHANK.Flatten.reconstruct(m_par, par_prior);
e_set = BASEforHANK.e_set;

# set some paths
@set! e_set.save_mode_file = paths["bld_example"] * "/HANK_mode.jld2";
@set! e_set.save_posterior_file = paths["bld_example"] * "/HANK_chain.jld2";
@set! e_set.mode_start_file = paths["src_example"] * "/Data/par_final_dict.txt";
@set! e_set.data_file = paths["src_example"] * "/Data/bbl_data_inequality.csv";

# fix seed for random number generation
BASEforHANK.Random.seed!(e_set.seed);

## ------------------------------------------------------------------------------------------
## Calculate Steady State and prepare linearization
## ------------------------------------------------------------------------------------------

# steady state at the prior mode
ss_full = call_find_steadystate(m_par);

# sparse DCT representation
sr_full = call_prepare_linearization(ss_full, m_par);

# save the steady state
jldsave(paths["bld_example"] * "/steadystate.jld2", true; sr_full);

# compute steady state moments
K = exp.(sr_full.XSS[sr_full.indexes.KSS]);
B = exp.(sr_full.XSS[sr_full.indexes.BSS]);
Bgov = exp.(sr_full.XSS[sr_full.indexes.BgovSS]);
Y = exp.(sr_full.XSS[sr_full.indexes.YSS]);
T10W = exp(sr_full.XSS[sr_full.indexes.TOP10WshareSS]);
G = exp.(sr_full.XSS[sr_full.indexes.GSS]);
distr_b = sum(sr_full.distrSS; dims = (2, 3))[:];
fr_borr = sum(distr_b[sr_full.n_par.grid_b .< 0]);

# Display steady state moments
@printf "\n"
pretty_table(
    [
        "Liquid to Illiquid Assets Ratio" B/K
        "Capital to Output Ratio" K / Y/4.0
        "Government Debt to Output Ratio" Bgov / Y/4.0
        "Government Spending to Output Ratio" G/Y
        "TOP 10 Wealth Share" T10W
        "Fraction of Borrower" fr_borr
    ];
    header = ["Variable", "Value"],
    title = "Steady State Moments",
    formatters = ft_printf("%.4f"),
)

## ------------------------------------------------------------------------------------------
## Linearize the full model, find sparse state-space representation
## ------------------------------------------------------------------------------------------

lr_full = linearize_full_model(sr_full, m_par);

# save the linearization
jldsave(paths["bld_example"] * "/linearresults.jld2", true; lr_full);

# sparse state-space representation
sr_reduc = model_reduction(sr_full, lr_full, m_par);
lr_reduc = update_model(sr_reduc, lr_full, m_par);

# save the reduction
jldsave(paths["bld_example"] * "/reduction.jld2", true; sr_reduc, lr_reduc);

# model timing
@printf "One model solution takes: \n"
@set! sr_reduc.n_par.verbose = false;
@btime lr_reduc = update_model(sr_reduc, lr_full, m_par);
@set! sr_reduc.n_par.verbose = true;

## ------------------------------------------------------------------------------------------
## Estimation
## ------------------------------------------------------------------------------------------

if e_set.estimate_model == true
    @printf "\n"
    @printf "Estimation...\n"

    # warning: estimation might take a long time!
    er_mode, posterior_mode, smoother_mode, sr_mode, lr_mode, m_par_mode =
        find_mode(sr_reduc, lr_reduc, m_par, e_set)

    # Only relevant output for later plotting will be saved.
    # If you require all smoother output including the variance estimates
    # over time, items 4 and 5, comment out the next line.
    # This increases the hard disk storage significantly.
    smoother_mode = (0.0, 0.0, smoother_mode[3], 0.0, 0.0, smoother_mode[6], 0.0)

    # Stores mode finding results in file e_set.save_mode_file
    jldsave(
        e_set.save_mode_file,
        true;
        posterior_mode,
        smoother_mode,
        sr_mode,
        lr_mode,
        er_mode,
        m_par_mode,
        e_set,
    )
    # !! warning: the provided mode file does not contain smoothed covars (smoother_mode[4] and [5])!!
    # @load BASEforHANK.e_set.save_mode_file posterior_mode sr_mode lr_mode er_mode m_par_mode smoother_mode

    sr_mc,
    lr_mc,
    er_mc,
    m_par_mc,
    draws_raw,
    posterior,
    accept_rate,
    par_final,
    hessian_sym,
    smoother_output = sample_posterior(sr_mode, lr_mode, er_mode, m_par_mode, e_set)

    # Only relevant output for later plotting will be saved.
    # If you want all smoother output including the variance estimates
    # over time, items 4 and 5, comment out the next line.
    # This increases the hard disk storage significantly.
    smoother_output = (0.0, 0.0, smoother_output[3], 0.0, 0.0, smoother_output[6], 0.0)

    # Stores mcmc results in file e_set.save_posterior_file
    jldsave(
        e_set.save_posterior_file,
        true;
        sr_mc,
        lr_mc,
        er_mc,
        m_par_mc,
        draws_raw,
        posterior,
        accept_rate,
        par_final,
        hessian_sym,
        smoother_output,
        e_set,
    )
    # !! The following file is not provided !!
    #      @load BASEforHANK.e_set.save_posterior_file sr_mc lr_mc er_mc  m_par_mc draws_raw posterior accept_rate par_final hessian_sym smoother_output e_set

    @printf "Estimation... Done. \n"
end

## ------------------------------------------------------------------------------------------
## Compute all IRFs, VDs, BCVDs, and historical decompositions
## ------------------------------------------------------------------------------------------

@printf "\n"
@printf "Compute IRFs, VDs, and BCVDs...\n"

# Get indices of the shocks
exovars = [getfield(sr_mc.indexes_r, shock_names[i]) for i = 1:length(shock_names)];

# Get standard deviations of the shocks
stds_mc = [getfield(m_par_mc, Symbol("σ_", i)) for i in shock_names];
stds_mode = [getfield(m_par_mode, Symbol("σ_", i)) for i in shock_names];

# Compute IRFs
IRFs_mc, _, IRFs_order = compute_irfs(
    exovars,
    lr_mc.State2Control,
    lr_mc.LOMstate,
    sr_mc.XSS,
    sr_mc.indexes_r;
    init_val = stds_mc,
);
IRFs_mode, _, _, = compute_irfs(
    exovars,
    lr_mode.State2Control,
    lr_mode.LOMstate,
    sr_mode.XSS,
    sr_mc.indexes_r;
    init_val = stds_mode,
);

# Compute variance decomposition of IRFs
VDs_mc = compute_vardecomp(IRFs_mc);
VDs_mode = compute_vardecomp(IRFs_mode);

# Compute business cycle frequency variance decomposition
VDbcs_mc, UnconditionalVar_mc =
    compute_vardecomp_bcfreq(exovars, stds_mc, lr_mc.State2Control, lr_mc.LOMstate);
VDbcs_mode, UnconditionalVar_mode =
    compute_vardecomp_bcfreq(exovars, stds_mode, lr_mode.State2Control, lr_mode.LOMstate);

# Compute historical decompositions
ShockContr, ShockContr_order = compute_hist_decomp(
    exovars,
    lr_mc.State2Control,
    lr_mc.LOMstate,
    smoother_output,
    sr_mc.indexes_r,
);

## ------------------------------------------------------------------------------------------
## Graphical outputs
## ------------------------------------------------------------------------------------------

@printf "\n"
@printf "Plotting...\n"

mkpath(paths["bld_example"] * "/IRFs");
plot_irfs(
    [
        (:Z, "TFP"),
        (:ZI, "Inv.-spec. tech."),
        (:μ, "Price markup"),
        (:μw, "Wage markup"),
        (:A, "Risk premium"),
        (:Rshock, "Mon. policy"),
        (:Gshock, "Structural deficit"),
        (:Tprogshock, "Tax progr."),
        (:Sshock, "Income risk"),
    ],
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    [(IRFs_mc, "HANK Mode"), (IRFs_mode, "HANK Posterior")],
    IRFs_order,
    sr_mc.indexes_r;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs",
    yscale = "standard",
    style_options = (lw = 2, color = [:blue, :red], linestyle = [:solid, :dash]),
)

mkpath(paths["bld_example"] * "/IRFs_cat");
plot_irfs_cat(
    Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    IRFs_mc,
    IRFs_order,
    sr_mc.indexes_r;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs_cat",
    yscale = "standard",
    style_options = (lw = 2, color = [:blue, :red], linestyle = [:solid, :dash]),
)

mkpath(paths["bld_example"] * "/VDs");
plot_vardecomp(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    [(VDs_mc, "HANK Mode"), (VDs_mode, "HANK Posterior")],
    IRFs_order,
    sr_mc.indexes_r;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDs",
)

mkpath(paths["bld_example"] * "/VDs_cat");
plot_vardecomp(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    [(VDs_mc, "HANK Mode"), (VDs_mode, "HANK Posterior")],
    IRFs_order,
    sr_mc.indexes_r;
    shock_categories = Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDs_cat",
)

mkpath(paths["bld_example"] * "/VDbcs");
plot_vardecomp_bcfreq(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    [(VDbcs_mc, "HANK Mode"), (VDbcs_mode, "HANK Posterior")],
    IRFs_order,
    sr_mc.indexes_r;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDbcs",
)

mkpath(paths["bld_example"] * "/VDbcs_cat");
plot_vardecomp_bcfreq(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    [(VDbcs_mc, "HANK Mode"), (VDbcs_mode, "HANK Posterior")],
    IRFs_order,
    sr_mc.indexes_r;
    shock_categories = Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDbcs_cat",
)

mkpath(paths["bld_example"] * "/HDs");
plot_hist_decomp(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    ShockContr,
    ShockContr_order,
    sr_mc.indexes_r;
    timeline = collect(1954.75:0.25:2019.75),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/HDs",
);

mkpath(paths["bld_example"] * "/HDs_cat");
plot_hist_decomp(
    [
        (:Ygrowth, "Output growth"),
        (:Cgrowth, "Consumption growth"),
        (:Igrowth, "Investment growth"),
        (:N, "Employment"),
        (:wgrowth, "Wage growth"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:σ, "Income risk"),
        (:Tprog, "Tax progressivity"),
        (:TOP10Wshare, "Top 10 wealth share"),
        (:TOP10Ishare, "Top 10 inc. share"),
    ],
    ShockContr,
    ShockContr_order,
    sr_mc.indexes_r;
    shock_categories = Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    timeline = collect(1954.75:0.25:2019.75),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/HDs_cat",
);

@printf "\n"
@printf "Done.\n"
