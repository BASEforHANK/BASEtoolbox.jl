"""
Mainboard for testing the results of the baseline example.
    - steady state of the full model
"""

debug = false;
include("preamble.jl");

## ------------------------------------------------------------------------------------------
## Calculate Steady State and prepare linearization
## ------------------------------------------------------------------------------------------

m_par = ModelParameters();

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

# testing the parameters
include("tests/tests_parameters.jl");

# steady state at m_par (the stable one – so that no carry over effects are present)
ss_full = call_find_steadystate(test_m_par);

# testing the steady state
include("tests/tests_steady_state.jl");

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
