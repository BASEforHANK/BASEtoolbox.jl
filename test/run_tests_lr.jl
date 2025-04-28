"""
Mainboard for testing the results of the baseline example.
    - linearization of the full model
"""

debug = false;
include("preamble.jl");

## ------------------------------------------------------------------------------------------
## Linearize the full model, find sparse state-space representation
## ------------------------------------------------------------------------------------------

lr_full = linearize_full_model(test_sr_full, test_m_par);

# testing the linearization
include("tests/tests_linearization.jl");

## ------------------------------------------------------------------------------------------
## Final check: did all tests pass?
## ------------------------------------------------------------------------------------------

if !all_passed
    @error "Some tests failed! Please check the output."
    exit(1) # Exit with a nonzero code to make GitHub Actions fail
end
