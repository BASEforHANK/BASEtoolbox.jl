"""
Mainboard for a one asset version of the baseline example, no estimation.
"""

using PrettyTables, Printf;

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
    "bld_example" => replace(@__DIR__, "examples" => "bld") * "_noestim",
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
## Initialize: set up model parameters
## ------------------------------------------------------------------------------------------

m_par = ModelParameters();
e_set = BASEforHANK.e_set;

@set! m_par.ξ = 4.0;
@set! m_par.γ = 2.0;
@set! m_par.β = 0.98255;
@set! m_par.λ = 0.065;
@set! m_par.ρ_h = 0.98;
@set! m_par.σ_h = 0.12;
@set! m_par.ι = 0.0625;
@set! m_par.ζ = 0.00022222222222222223;
@set! m_par.α = 0.318;
@set! m_par.δ_0 = 0.021500000000000002;
@set! m_par.δ_s = 0.7055720197078786;
@set! m_par.ϕ = 1.9409223183717077;
@set! m_par.μ = 1.1;
@set! m_par.κ = 0.1456082664986374;
@set! m_par.μw = 1.1;
@set! m_par.κw = 0.23931075416274708;
@set! m_par.Tlev = 1.0 + (1 - 0.8225);
@set! m_par.Tprog = 1.0 + 0.1022;
@set! m_par.RRB = 1.0;
@set! m_par.Rbar = 0.021778180864641117;
@set! m_par.ωΠ = 0.2;
@set! m_par.ιΠ = 0.016;
@set! m_par.shiftΠ = 0.7002848330469671;
@set! m_par.ρ_A = 0.9724112284399131;
@set! m_par.σ_A = 0.0015812471705012755;
@set! m_par.ρ_Z = 0.9978155269262137;
@set! m_par.σ_Z = 0.00600947811158941;
@set! m_par.ρ_ZI = 0.7637111671257767;
@set! m_par.σ_ZI = 0.0721141538701523;
@set! m_par.ρ_μ = 0.903740078830077;
@set! m_par.σ_μ = 0.01350860622318172;
@set! m_par.ρ_μw = 0.9057892147641305;
@set! m_par.σ_μw = 0.035058308969408175;
@set! m_par.ρ_s = 0.544722245741144;
@set! m_par.σ_Sshock = 0.6918558038597916;
@set! m_par.Σ_n = 28.879770107327673;
@set! m_par.ρ_R = 0.8030565250630299;
@set! m_par.σ_Rshock = 0.002306627917745612;
@set! m_par.θ_π = 2.0780841671981856;
@set! m_par.θ_Y = 0.21872568927661648;
@set! m_par.γ_B = 0.020131162775595176;
@set! m_par.γ_π = -2.1737350397931947;
@set! m_par.γ_Y = -0.4363130165391906;
@set! m_par.ρ_Gshock = 0.9682224473297878;
@set! m_par.σ_Gshock = 0.003761816459554433;
@set! m_par.ρ_τ = 0.4926482696848203;
@set! m_par.γ_Bτ = 3.293063617271948;
@set! m_par.γ_Yτ = -0.9207283604196101;
@set! m_par.ρ_P = 0.9194235885358465;
@set! m_par.σ_Tprogshock = 0.06865440038519788;
@set! m_par.γ_BP = 0.0;
@set! m_par.γ_YP = 0.0;
@set! m_par.γ_WP = 0.0;
@set! m_par.ρ_Rshock = 1.0e-8;
@set! m_par.ρ_Tprogshock = 1.0e-8;
@set! m_par.ρ_Sshock = 1.0e-8;

## ------------------------------------------------------------------------------------------
## Calculate Steady State and prepare linearization
## ------------------------------------------------------------------------------------------

# steady state at m_par
ss_full = call_find_steadystate(m_par);

# sparse DCT representation
sr_full = call_prepare_linearization(ss_full, m_par);

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

## ------------------------------------------------------------------------------------------
## Compute all IRFs, VDs, and BCVDs
## ------------------------------------------------------------------------------------------

@printf "\n"
@printf "Compute IRFs, VDs, and BCVDs...\n"

# Get indices of the shocks
exovars = [getfield(sr_full.indexes, shock_names[i]) for i = 1:length(shock_names)];

# Get standard deviations of the shocks
stds = [getfield(sr_full.m_par, Symbol("σ_", i)) for i in shock_names];

# Compute IRFs
IRFs, _, IRFs_order = compute_irfs(
    exovars,
    lr_full.State2Control,
    lr_full.LOMstate,
    sr_full.XSS,
    sr_full.indexes;
    init_val = stds,
);

# Compute variance decomposition of IRFs
VDs = compute_vardecomp(IRFs);

# Compute business cycle frequency variance decomposition
VDbcs, UnconditionalVar =
    compute_vardecomp_bcfreq(exovars, stds, lr_full.State2Control, lr_full.LOMstate);

# Compute distributional IRFs
_, _, _, IRFs_dist = compute_irfs(
    exovars,
    lr_full.State2Control,
    lr_full.LOMstate,
    sr_full.XSS,
    sr_full.indexes;
    init_val = stds,
    distribution = true,
    comp_ids = sr_full.compressionIndexes,
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
    [(IRFs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs",
    yscale = "standard",
    style_options = (lw = 2, color = [:blue, :red], linestyle = [:solid, :dash]),
);

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
    IRFs,
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs_cat",
    yscale = "standard",
    style_options = (lw = 2, color = [:blue, :red], linestyle = [:solid, :dash]),
);

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
    [(VDs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDs",
);

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
    [(VDs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    shock_categories = Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDs_cat",
);

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
    [(VDbcs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDbcs",
);

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
    [(VDbcs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    shock_categories = Dict(
        ("Monetary", "mon") => [:Rshock, :A],
        ("Fiscal", "fis") => [:Gshock, :Tprogshock],
        ("Productivity", "pro") => [:Z, :ZI, :μ, :μw],
    ),
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDbcs_cat",
);

mkpath(paths["bld_example"] * "/IRFs_dist");
plot_distributional_irfs(
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
        ("Wb", "Marginal Value of Bonds"),
        ("Wk", "Marginal Value of Capital"),
        ("distr_b", "Bonds"),
        ("distr_k", "Capital"),
    ],
    IRFs_dist,
    IRFs_order;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs_dist",
);

@printf "\n"
@printf "Done.\n"
