# This file defines the sets of aggregate shocks, states (inluding shocks), and controls
# The code checks that the number of equations in the aggregate model is equal to the number
# of aggregate variables excluding the distributional summary statistics. The latter are not
# contained in the aggregate model code as they are parameter free but change whenever the
# distribution changes and do not show up in any aggregate model equation.

n_rep = 1 # Number of n_rep of some model equations (e.g. countries, industries)

# List of aggregate shocks, without duplication (e.g. across countries or industries)
shock_names = [:Z, :ZI, :μ, :μw, :A, :Rshock, :Gshock, :Tprogshock, :Sshock]

# List of aggregate states, without duplication of names (e.g. across countries or industries)
# Duplicated names are created below
state_names = [
    "A",
    "Z",
    "ZI",
    "RB",
    "μ",
    "μw",
    "σ",
    "Ylag",
    "Bgovlag",
    "Tlag",
    "Ilag",
    "wFlag",
    "qlag",
    "Clag",
    "Tbarlag",
    "Tproglag",
    "qΠlag",
    "Gshock",
    "Tprogshock",
    "Rshock",
    "Sshock",
    "K",
    "B",
]

# List of (the subset) of aggregate states, that need to be duplicated (e.g. across countries or industries)
dup_state_names = Vector{String}(undef, 0)
#[
#     "A",
#]

# List cross-sectional controls / distributional summary variables (no equations in aggregate model expected)
# if these need to be duplicated, do by hand !
distr_names = ["GiniC", "GiniW", "TOP10Ishare", "TOP10Inetshare", "TOP10Wshare", "sdlogy"]

control_names = [
    "RK",
    "wF",
    "π",
    "πw",
    "Y",
    "C",
    "q",
    "N",
    "mc",
    "mcw",
    "u",
    "qΠ",
    "Π_F",
    "RL",
    "RD",
    "RRL",
    "RRD",
    "Bgov",
    "Hprog",
    "Htilde",
    "Tbar",
    "T",
    "I",
    "BY",
    "TY",
    "wH",
    "G",
    "Tlev",
    "Tprog",
    "Tc",
    "Ygrowth",
    "Bgovgrowth",
    "Igrowth",
    "wgrowth",
    "Cgrowth",
    "Tgrowth",
    "LP",
    "LPXA",
    "Π_U",
    "Π_E",
    "TotalAssets",
]

# List of (the subset) of aggregate states, that need to be duplicated (e.g. across countries or industries)
dup_control_names = Vector{String}(undef, 0)#
#["",
#     "LP",
#]

# Delete duplicated state and shock names and create duplicated variable names of state and control names

state_names = setdiff(state_names, dup_state_names)
for j = 1:n_rep # Create names of duplicated state variables (e.g. across countries or industries)
    if j == 1
        aux = ""
    else
        aux = string(j)
    end
    append!(state_names, dup_state_names .* aux)
end

shock_states = string.(shock_names)
dup_shock_states = intersect(string.(shock_names), dup_state_names)
shock_states = setdiff(shock_states, dup_shock_states)
for j = 1:n_rep # Create names of duplicated shock variables (e.g. across countries or industries)
    if j == 1
        aux = ""
    else
        aux = string(j)
    end
    append!(shock_states, dup_shock_states .* aux)
end
shock_names = Symbol.(shock_states)

control_names = setdiff(control_names, dup_control_names)
for j = 1:n_rep # Create names of duplicated control variables (e.g. across countries or industries)
    if j == 1
        aux = ""
    else
        aux = string(j)
    end
    append!(control_names, dup_control_names .* aux)
end

# All controls in one array
control_names = [distr_names; control_names]
# All names in one array
aggr_names = [state_names; control_names]

args_hh_prob_names = [
    "wH",
    "N",
    "Hprog",
    "q",
    "RRL",
    "RRD",
    "RK",
    "Tlev",
    "Tprog",
    "Tbar",
    "Tc",
    "Π_E",
    "Π_U",
    "Htilde",
    "σ",
]

# ascii names used for cases where unicode doesn't work, e.g., file saves
unicode2ascii(x) = replace.(
    replace.(
        replace.(replace.(replace.(x, "τ" => "tau"), "σ" => "sigma"), "π" => "pi"),
        "μ" => "mu",
    ),
    "ρ" => "rho",
)

state_names_ascii = unicode2ascii(state_names)
control_names_ascii = unicode2ascii(control_names)
aggr_names_ascii = [state_names_ascii; control_names_ascii]
