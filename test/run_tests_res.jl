"""
Mainboard for testing the results of the baseline example.
    - results of the full model
"""

debug = false;
include("preamble.jl");

## ------------------------------------------------------------------------------------------
## Graphical outputs
## ------------------------------------------------------------------------------------------

shock_names = [:Z, :ZI, :μ, :μw, :A, :Rshock, :Gshock, :Tprogshock, :Sshock];

# Get indices of the shocks
exovars = [getfield(test_sr_full.indexes, shock_names[i]) for i = 1:length(shock_names)];

# Get standard deviations of the shocks
stds = [getfield(test_sr_full.m_par, Symbol("σ_", i)) for i in shock_names];

# Compute IRFs
IRFs, _, IRFs_order = compute_irfs(
    exovars,
    test_lr_full.State2Control,
    test_lr_full.LOMstate,
    test_sr_full.XSS,
    test_sr_full.indexes;
    init_val = stds,
);

# Compute variance decomposition of IRFs
VDs = compute_vardecomp(IRFs);

# Compute business cycle frequency variance decomposition
VDbcs, UnconditionalVar = compute_vardecomp_bcfreq(
    exovars,
    stds,
    test_lr_full.State2Control,
    test_lr_full.LOMstate,
);

results = Dict("IRFs" => IRFs[:, 1:100, :], "VDs" => VDs[:, 1:100, :], "VDbcs" => VDbcs);

# testing the results
include("tests/tests_results.jl")

## ------------------------------------------------------------------------------------------
## Final check: did all tests pass?
## ------------------------------------------------------------------------------------------

if !all_passed
    @error "Some tests failed! Please check the output."
    exit(1) # Exit with a nonzero code to make GitHub Actions fail
end
