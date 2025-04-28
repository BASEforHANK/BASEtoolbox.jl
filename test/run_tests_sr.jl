"""
Mainboard for testing the results of the baseline example.
    - prepare linearization of the full model
"""

debug = false;
include("preamble.jl");

## ------------------------------------------------------------------------------------------
## Calculate Steady State and prepare linearization
## ------------------------------------------------------------------------------------------

# sparse DCT representation (the stable one – so that no carry over effects are present)
sr_full = call_prepare_linearization(test_ss_full, test_m_par);

# testing the reduction
include("tests/tests_reduction.jl");

## ------------------------------------------------------------------------------------------
## Final check: did all tests pass?
## ------------------------------------------------------------------------------------------

if !all_passed
    @error "Some tests failed! Please check the output."
    exit(1) # Exit with a nonzero code to make GitHub Actions fail
end
