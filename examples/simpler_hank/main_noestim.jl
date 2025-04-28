"""
Mainboard for the example "simpler HANK model" of the BASEforHANK package, no estimation.
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

## ------------------------------------------------------------------------------------------
## Graphical outputs
## ------------------------------------------------------------------------------------------

@printf "\n"
@printf "Plotting...\n"

mkpath(paths["bld_example"] * "/IRFs");
plot_irfs(
    [(:Z, "TFP"), (:Rshock, "Mon. policy")],
    [
        (:Y, "Output"),
        (:C, "Consumption"),
        (:I, "Investment"),
        (:N, "Employment"),
        (:wF, "Wage (firms)"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:B, "Government debt"),
        (:Tprog, "Tax progressivity"),
        (:Tlev, "Tax level"),
        (:G, "Government spending"),
    ],
    [(IRFs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/IRFs",
);

mkpath(paths["bld_example"] * "/VDs");
plot_vardecomp(
    [
        (:Y, "Output"),
        (:C, "Consumption"),
        (:I, "Investment"),
        (:N, "Employment"),
        (:wF, "Wage (firms)"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:B, "Government debt"),
        (:Tprog, "Tax progressivity"),
        (:Tlev, "Tax level"),
        (:G, "Government spending"),
    ],
    [(VDs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDs",
);

mkpath(paths["bld_example"] * "/VDbcs");
plot_vardecomp_bcfreq(
    [
        (:Y, "Output"),
        (:C, "Consumption"),
        (:I, "Investment"),
        (:N, "Employment"),
        (:wF, "Wage (firms)"),
        (:RB, "Nominal rate"),
        (:π, "Inflation"),
        (:B, "Government debt"),
        (:Tprog, "Tax progressivity"),
        (:Tlev, "Tax level"),
        (:G, "Government spending"),
    ],
    [(VDbcs, "Baseline")],
    IRFs_order,
    sr_full.indexes;
    show_fig = false,
    save_fig = true,
    path = paths["bld_example"] * "/VDbcs",
);

@printf "\n"
@printf "Done.\n"
