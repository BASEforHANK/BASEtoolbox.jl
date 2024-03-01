#------------------------------------------------------------------------------
# Header: load module
#------------------------------------------------------------------------------
# ATTENTION: make sure that your present working directory pwd() is set to the folder
# containing script.jl and BASEforHANK.jl. Otherwise adjust the load path.
cd("./src")
# push!(LOAD_PATH, pwd())
# pre-process user inputs for model setup
include("Preprocessor/PreprocessInputs.jl")
include("BASEforHANK.jl")
using .BASEforHANK
using BenchmarkTools
# set BLAS threads to the number of Julia threads.
# prevents BLAS from grabbing all threads on a machine
BASEforHANK.LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

#------------------------------------------------------------------------------
# initialize parameters to priors to select coefficients of DCTs of Vm, Vk]
# that are retained 
#------------------------------------------------------------------------------
m_par = ModelParameters()
priors = collect(metaflatten(m_par, prior)) # model parameters

par_prior = mode.(priors)
m_par = BASEforHANK.Flatten.reconstruct(m_par, par_prior)
e_set = BASEforHANK.e_set;
# alternatively, load estimated parameters by running, e.g.,
# @load BASEforHANK.e_set.save_posterior_file par_final e_set
# m_par = BASEforHANK.Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(e_set.meas_error_input)])

# Fix seed for random number generation
BASEforHANK.Random.seed!(e_set.seed)

################################################################################
# Comment in the following block to be able to go straight to plotting (comment out lines 40-53)
################################################################################
@load "Output/Saves/steadystate.jld2" sr_full
@load "Output/Saves/linearresults.jld2" lr_full
@load "Output/Saves/reduction.jld2" sr_reduc lr_reduc
# @load BASEforHANK.e_set.save_posterior_file sr_mc lr_mc er_mc m_par_mc smoother_output
# @set! e_set.estimate_model = false 

# Calculate Steady State at prior mode 
println("Calculating the steady state")
ss_full = call_find_steadystate(m_par)
# Find sparse DCT representation
println("preparing the linearization")
sr_full = call_prepare_linearization(ss_full, m_par)
# COMPACT call of both of the above:
# sr_full = compute_steadystate(m_par)

jldsave("Output/Saves/steadystate.jld2", true; sr_full) # true enables compression
# @load "Output/Saves/steadystate.jld2" sr_full

#------------------------------------------------------------------------------
# compute and display steady-state moments
#------------------------------------------------------------------------------
K = exp.(sr_full.XSS[sr_full.indexes.KSS])
B = exp.(sr_full.XSS[sr_full.indexes.BSS])
Bgov = exp.(sr_full.XSS[sr_full.indexes.BgovSS])
Y = exp.(sr_full.XSS[sr_full.indexes.YSS])
T10W = exp(sr_full.XSS[sr_full.indexes.TOP10WshareSS])
G = exp.(sr_full.XSS[sr_full.indexes.GSS])
distr_m = sum(sr_full.distrSS, dims = (2, 3))[:]
fr_borr = sum(distr_m[sr_full.n_par.grid_m.<0])

println("Steady State Moments:")
println("Liquid to Illiquid Assets Ratio:", B / K)
println("Capital to Output Ratio:", K / Y / 4.0)
println("Government Debt to Output Ratio:", Bgov / Y / 4.0)
println("Government spending to Output Ratio:", G / Y)
println("TOP 10 Wealth Share:", T10W)
println("Fraction of Borrower:", fr_borr)

# linearize the full model
lr_full = linearize_full_model(sr_full, m_par)
jldsave("Output/Saves/linearresults.jld2", true; lr_full)
# @load "Output/Saves/linearresults.jld2" lr_full

# Find sparse state-space representation
sr_reduc = model_reduction(sr_full, lr_full, m_par);
lr_reduc = update_model(sr_reduc, lr_full, m_par)
jldsave("Output/Saves/reduction.jld2", true; sr_reduc, lr_reduc)
# @load "Output/Saves/reduction.jld2" sr_reduc lr_reduc

# model timing
println("One model solution takes")
@set! sr_reduc.n_par.verbose = false
@btime lr_reduc = update_model(sr_reduc, lr_full, m_par)
@set! sr_reduc.n_par.verbose = true;

if e_set.estimate_model == true

    # warning: estimation might take a long time!
    er_mode, posterior_mode, smoother_mode, sr_mode, lr_mode, m_par_mode =
        find_mode(sr_reduc, lr_reduc, m_par)
    
    # Only relevant output for later plotting will be saved.
    # If you require all smoother output including the variance estimates 
    # over time, items 4 and 5, comment out the next line.  
    # This increases the hard disk storage significantly.
    smoother_mode = (0.0,0.0,smoother_mode[3],0.0,0.0, smoother_mode[6],0.0)

    # Stores mode finding results in file e_set.save_mode_file 
    jldsave(
        BASEforHANK.e_set.save_mode_file,
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
    smoother_output = sample_posterior(sr_mode, lr_mode, er_mode, m_par_mode)
    
    # Only relevant output for later plotting will be saved.
    # If you want all smoother output including the variance estimates 
    # over time, items 4 and 5, comment out the next line.  
    # This increases the hard disk storage significantly.
    smoother_output = (0.0,0.0,smoother_output[3],0.0,0.0, smoother_output[6],0.0)
    
    # Stores mcmc results in file e_set.save_posterior_file 
    jldsave(
        BASEforHANK.e_set.save_posterior_file,
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

end


##############################################################################################
# Graphical Model Output
###############################################################################################
# variables to be plotted
select_variables = [
    :Ygrowth,
    :Cgrowth,
    :Igrowth,
    :N,
    :wgrowth,
    :RB,
    :π,
    :σ,
    :τprog,
    :TOP10Wshare,
    :TOP10Ishare,
]

# models to be plotted
number_models = 2
model_names = Array{String}(undef, 1, number_models)
model_names[1] = "HANK Mode"
model_names[2] = "HANK Posterior"

# enter here the models, as tupel of tupels (sr, lr, e_set, m_par), to be compared
models_tupel = ((sr_mc, lr_mc, e_set, m_par_mc), (sr_mode, lr_mode, e_set, m_par_mode))

timeline = collect(1954.75:0.25:2019.75)
select_vd_horizons = [4 16 100] # horizons for variance decompositions
recessions_vec = [
    1957.5,
    1958.25,
    1960.25,
    1961.0,
    1969.75,
    1970.75,
    1973.75,
    1975.0,
    1980.0,
    1980.5,
    1981.5,
    1982.75,
    1990.5,
    1991.0,
    2001.0,
    2001.75,
    2007.75,
    2009.25,
] # US recession dates for plotting

# "nice" names for labels
nice_var_names = [
    "Output growth",
    "Consumption growth",
    "Investment growth",
    "Employment",
    "Wage growth",
    "Nominal rate",
    "Inflation",
    "Income risk",
    "Tax progressivity",
    "Top 10 wealth share",
    "Top 10 inc. share",
]
nice_s_names = [
    "TFP",
    "Inv.-spec. tech.",
    "Price markup",
    "Wage markup",
    "Risk premium",
    "Mon. policy",
    "Structural deficit",
    "Tax progr.",
    "Income risk",
]

# compute IRFs for all models in tupel, all variables in select_variables
IRFs, VDs, SHOCKs, VD_bc_s = compute_irfs_vardecomp(models_tupel, select_variables)

# display IRFs and export as pdf
IRFs_plot = plot_irfs(
    IRFs,
    SHOCKs,
    select_variables,
    nice_var_names,
    nice_s_names,
    40,
    model_names,
    4;
    savepdf = true,
    disp_switch = true
)

# export Variance Decompositions as DataFrames and Plot using VegaLite
DF_V_Decomp, DF_V_Decomp_bc = plot_vardecomp(
    VDs,
    VD_bc_s,
    select_vd_horizons,
    model_names,
    SHOCKs,
    select_variables;
    savepdf = true,
    suffix = "_nolegend",
    legend_switch = true,
    disp_switch = true
)

# produce historical contributions as Array and Data Frame and plot p
Historical_contrib_HA, DF_H_Decomp_HA, HD_plot_HA = compute_hist_decomp(
    sr_mc,
    lr_mc,
    e_set,
    m_par_mc,
    smoother_output,
    select_variables,
    timeline;
    savepdf = true,
    prefix = "HA_",
)
