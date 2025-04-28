"""
    @generate_equations()

Write out the expansions around steady state for all variables in `aggr_names`, i.e.
generate code that reads aggregate states/controls from steady state deviations.

Equations take the form of (with variable `r` as example):

  - `r       = exp.(XSS[indexes.rSS] .+ X[indexes.r])`
  - `rPrime  = exp.(XSS[indexes.rSS] .+ XPrime[indexes.r])`

# Requires

(module) global `aggr_names`
"""
macro generate_equations()
    ex = quote end
    for j in aggr_names
        i = Symbol(j)
        varnamePrime = Symbol(i, "Prime")
        varnameSS = Symbol(i, "SS")
        ex_aux = quote
            $i = exp.(XSS[indexes.$varnameSS] .+ X[indexes.$i])
            $varnamePrime = exp.(XSS[indexes.$varnameSS] .+ XPrime[indexes.$i])
        end
        append!(ex.args, ex_aux.args)
    end
    return esc(ex)
end

"""
    @write_args_hh_prob()

Write all variables defined in the global string vector `args_hh_prob_names` into a vector
`args_hh_prob` based on the variables with the according names in local scope. The type of
the vector is inferred from the first variable in `args_hh_prob_names`.

# Requires

(module) global `args_hh_prob_names`
"""
macro write_args_hh_prob()
    varname = Symbol(args_hh_prob_names[1])
    ex = quote
        args_hh_prob = Vector{typeof($varname)}(undef, length(args_hh_prob_names))
    end
    for (i, j) in enumerate(args_hh_prob_names)
        varname = Symbol(j)
        ex_aux = quote
            args_hh_prob[$i] = $varname
        end
        append!(ex.args, ex_aux.args)
    end
    return esc(ex)
end

"""
    @write_args_hh_prob_ss()

See [`@write_args_hh_prob`](@ref), however, for the case where the variables are called with
a suffix "SS".

# Requires

(module) global `args_hh_prob_names`
"""

macro write_args_hh_prob_ss()
    varname = Symbol(args_hh_prob_names[1], "SS")
    ex = quote
        args_hh_prob = Vector{typeof($varname)}(undef, length(args_hh_prob_names))
    end
    for (i, j) in enumerate(args_hh_prob_names)
        varname = Symbol(j, "SS")
        ex_aux = quote
            args_hh_prob[$i] = $varname
        end
        append!(ex.args, ex_aux.args)
    end
    return esc(ex)
end

"""
    @read_args_hh_prob()

Read all variables defined in the global string vector `args_hh_prob_names` into local scope
based on the variables with the according names in the vector `args_hh_prob`.

# Requires

(module) global `args_hh_prob_names`, `args_hh_prob`
"""
macro read_args_hh_prob()
    ex = quote end
    for (i, j) in enumerate(args_hh_prob_names)
        varname = Symbol(j)
        ex_aux = quote
            $varname = args_hh_prob[$i]
        end
        append!(ex.args, ex_aux.args)
    end
    return esc(ex)
end

"""
    @read_args_hh_prob_ss()

See [`@read_args_hh_prob`](@ref), however, for the case where the variables are called with
a suffix "SS".

# Requires

(module) global `args_hh_prob_names`, `args_hh_prob`
"""

macro read_args_hh_prob_ss()
    ex = quote end
    for (i, j) in enumerate(args_hh_prob_names)
        varname = Symbol(j, "SS")
        ex_aux = quote
            $varname = args_hh_prob[$i]
        end
        append!(ex.args, ex_aux.args)
    end
    return esc(ex)
end

"""
    @writeXSS()

Write all single steady state variables into vectors XSS / XSSaggr.

# Requires

(module) globals `state_names`, `control_names`, `aggr_names`
"""
macro writeXSS()
    ex = quote
        XSS = [distr_bSS[:]; distr_kSS[:]; distr_hSS[:]; distrSS[:]]
    end
    for j in state_names
        varnameSS = Symbol(j, "SS")
        ex_aux = quote
            append!(XSS, log($varnameSS))
        end
        append!(ex.args, ex_aux.args)
    end

    ex_aux = quote
        append!(XSS, [WbSS[:]; WkSS[:]])
    end
    append!(ex.args, ex_aux.args)
    for j in control_names
        varnameSS = Symbol(j, "SS")
        ex_aux = quote
            append!(XSS, log($varnameSS))
        end
        append!(ex.args, ex_aux.args)
    end
    ex_aux = quote
        XSSaggr = [0.0]
    end
    append!(ex.args, ex_aux.args)
    for j in aggr_names
        varnameSS = Symbol(j, "SS")
        ex_aux = quote
            append!(XSSaggr, log($varnameSS))
        end
        append!(ex.args, ex_aux.args)
    end
    ex_aux = quote
        deleteat!(XSSaggr, 1)
    end
    append!(ex.args, ex_aux.args)

    return esc(ex)
end

"""
    @make_fnaggr(fn_name)

Create function `fn_name` that returns an instance of `IndexStructAggr` (created by
[`@make_struct_aggr`](@ref)), mapping aggregate states and controls to values `1` to
`length(aggr_names)` (both steady state and deviation from it).

# Requires

(module) global `aggr_names`
"""
macro make_fnaggr(fn_name)
    n_states = length(aggr_names)
    fieldsSS_states = [:($i) for i = 1:n_states]
    fields_states = [:($i) for i = 1:n_states]
    esc(quote
        function $(fn_name)(n_par)
            indexes = IndexStructAggr($(fieldsSS_states...), $(fields_states...))
            return indexes
        end
    end)
end

"""
    @make_fn(fn_name)

Create function `fn_name` that returns an instance of `IndexStruct` (created by
[`@make_struct`](@ref)), mapping states and controls to indexes inferred from numerical
parameters and compression indexes.

# Requires

(module) global `state_names`, `control_names`
"""
macro make_fn(fn_name)
    n_states = length(state_names)
    n_controls = length(control_names)
    fieldsSS_states = [:((tNo + tNo2) + $i) for i = 1:n_states]
    fields_states = [:(tNo + tNo4 - 3 + $i) for i = 1:n_states]
    fieldsSS_controls = [:(tNo + 3 * tNo2 + $i) for i in n_states .+ (1:n_controls)]
    fields_controls = [:(tNo + tNo3 + tNo4 - 3 + $i) for i in n_states .+ (1:n_controls)]
    esc(
        quote
            function $(fn_name)(
                n_par,
                compressionIndexesWb,
                compressionIndexesWk,
                compressionIndexesD,
            )
                tNo = n_par.nb + n_par.nk + n_par.nh
                tNo2 = n_par.nb * n_par.nk * n_par.nh
                tNo3 = length(compressionIndexesWb) + length(compressionIndexesWk)
                tNo4 = length(compressionIndexesD)
                indexes = IndexStruct(
                    1:(n_par.nb), # distr_bSS
                    (n_par.nb + 1):(n_par.nb + n_par.nk), # distr_kSS
                    (n_par.nb + n_par.nk + 1):(tNo), # distr_hSS
                    (tNo + 1):(tNo + tNo2), # COPSS
                    $(fieldsSS_states...),
                    ((tNo + tNo2) + $(n_states) + 1):((tNo + tNo2) + tNo2 + $(n_states)), # WbSS
                    ((tNo + tNo2) + tNo2 + $(n_states) + 1):((tNo + tNo2) + 2 * tNo2 + $(n_states)), # WkSS
                    $(fieldsSS_controls...),
                    1:(n_par.nb - 1), # distr_b
                    (n_par.nb):(n_par.nb + n_par.nk - 2), # distr_k
                    (n_par.nb + n_par.nk - 1):(tNo - 3), # distr_h
                    (tNo - 2):(tNo + tNo4 - 3), # COP
                    $(fields_states...),
                    (tNo + tNo4 + $(n_states) - 2):(tNo + tNo4 + length(
                        compressionIndexesWb,
                    ) + $(n_states) - 3), # Wb
                    (tNo + tNo4 + length(
                        compressionIndexesWb,
                    ) + $(n_states) - 2):(tNo + tNo4 + length(
                        compressionIndexesWb,
                    ) + length(compressionIndexesWk) + $(n_states) - 3), # Wk
                    $(fields_controls...),
                )
                return indexes
            end
        end,
    )
end
