using Test, JLD2

all_passed = true;

# Switch of all reporting
if !debug
    oldstd = stdout
    redirect_stdout(devnull)
end;

root_dir = replace(Base.current_project(), "Project.toml" => "");
cd(root_dir);

# set up paths for the project
paths = Dict(
    "root" => root_dir,
    "src" => joinpath(root_dir, "src"),
    "bld" => joinpath(root_dir, "test/bld"),
    "src_example" => joinpath(root_dir, "examples/baseline"),
    "bld_example" => joinpath(root_dir, "test/bld/baseline"),
);

# create bld directory for the current example
mkpath(paths["bld_example"]);

# pre-process user inputs for model setup
include(paths["src"] * "/Preprocessor/PreprocessInputs.jl");
include(paths["src"] * "/BASEforHANK.jl");
using .BASEforHANK;

# set BLAS threads to the number of Julia threads, prevents grabbing all
BASEforHANK.LinearAlgebra.BLAS.set_num_threads(Threads.nthreads());

include("tests/helpers.jl");

## ------------------------------------------------------------------------------------------
## Load all test results
## ------------------------------------------------------------------------------------------

# Load the stored file(s)
test_m_par = load("test/results/m_par.jld2", "m_par");
test_ss_full = load("test/results/ss_full.jld2", "ss_full");
test_sr_full = load("test/results/sr_full.jld2", "sr_full");
test_lr_full = load("test/results/lr_full.jld2", "lr_full");
test_results = load("test/results/results.jld2");
