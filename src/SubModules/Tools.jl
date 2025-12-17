
module Tools

using ..Types

using LinearAlgebra,
    CategoricalArrays, Roots, ForwardDiff, Distributions, Printf, Dates, Logging
using SpecialFunctions: erf
using FFTW: dct, ifft

export Brent,
    broyden,
    CustomBrent,
    centralderiv,
    mydctmx,
    uncompress,
    compress,
    select_ind,
    Fastroot,
    my_integrate,
    myinterpolate3,
    mylinearinterpolate,
    mylinearinterpolate2,
    mylinearinterpolate3,
    mylinearinterpolate_mult2,
    mylinearinterpolate_mult3,
    mylinearinterpolate!,
    mylinearinterpolate_mult2!,
    mylinearinterpolate_mult2!,
    mylinearinterpolate_mult3!,
    locate,
    Tauchen,
    ExTransition,
    cdf_to_pdf,
    pdf_to_cdf,
    real_schur,
    dualpart,
    realpart,
    distrSummaries,
    gini,
    struc_to_vec,
    vec_to_struct,
    read_mem_linux,
    timer_help,
    stepwiseRSKron,
    stepwiseLSKron,
    doublingGxx,
    doublingSimple

include("Tools/BrentsMethod.jl")
include("Tools/LinInterpols.jl")
include("Tools/LocateFcn.jl")
include("Tools/GCintegration.jl")
include("Tools/DCT.jl")
include("Tools/FastRoot.jl")
include("Tools/MarkovChain.jl")
include("Tools/SchurUtils.jl")
include("Tools/DualUtils.jl")
include("Tools/Pdf2cdf.jl")
include("Tools/Broyden.jl")
include("Tools/CentralDerivatives.jl")
include("Tools/StructToVec.jl")
include("Tools/SystemUtils.jl")
include("Tools/fastgensylv.jl")
include("Tools/stepwiseKronecker2.jl")

const ORIG_STDOUT = stdout
const ORIG_STDERR = stderr

"""
    unmute_println(xs...)

Print the given arguments to the original stdout (Tools.ORIG_STDOUT), bypassing any
redirection of Base.stdout that may be in effect.

Arguments

  - `xs...`: Zero or more values to print; forwarded to `println`.

Returns

  - `nothing`. Side effect: writes a newline-terminated line to `Tools.ORIG_STDOUT`.

Notes

  - This is useful when output has been redirected (e.g., via `redirect_stdout`) but you
    still want to emit a message to the original terminal or log sink that was captured in
    `Tools.ORIG_STDOUT`.
  - Does not change the current stdout; only writes directly to `Tools.ORIG_STDOUT`.
"""
function unmute_println(xs...)
    println(Tools.ORIG_STDOUT, xs...)
end

"""
    unmute_printf(fmt::AbstractString, args...)

Write a formatted string to the original stdout (Tools.ORIG_STDOUT) using
`Printf`-style formatting, bypassing any redirection of Base.stdout.

Arguments

  - `fmt::AbstractString`: A `Printf`-style format string.
  - `args...`: Values to be formatted according to `fmt`.

Returns

  - `nothing`. Side effect: writes the formatted output (without adding an extra newline) to
    `Tools.ORIG_STDOUT`.

Notes

  - Equivalent in purpose to `@printf` or `print` but ensures the output goes to the
    preserved original stdout regardless of current redirection.
  - Formatting is performed via `Printf.Format(fmt)` semantics; invalid format strings or
    argument mismatches propagate the usual `Printf` errors.
"""
function unmute_printf(fmt::AbstractString, args...)
    print(Tools.ORIG_STDOUT, Printf.format(Printf.Format(fmt), args...))
end

"""
    quiet_call(f, args...; kwargs...)

Run `f(args...; kwargs...)` while suppressing all text written to `Base.stdout` and
`Base.stderr` by redirecting them to the system null device. The function returns the value
produced by `f` (or propagates any exception thrown by `f`).

Arguments

  - `f`: A callable to invoke.
  - `args...`: Positional arguments forwarded to `f`.
  - `kwargs...`: Keyword arguments forwarded to `f`.

Returns

  - Whatever `f` returns.

Notes

  - Only suppresses writes to `Base.stdout` and `Base.stderr` (i.e., low-level printing). It
    does not change the active logger; logging macros that use the global logger may still
    emit output unless the logger itself is changed.
  - Implementation uses the OS null device path (`"/dev/null"`), so behavior is
    platform-dependent (on Windows, use `"NUL"` instead).
  - Any exceptions raised by `f` propagate out of `quiet_call`; the redirections are
    properly restored when the call exits.
"""
function quiet_call(f, args...; kwargs...)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            return f(args...; kwargs...)
        end
    end
end

"""
    @silent expr

Evaluate `expr` while silencing both logging macros (by using a `Logging.NullLogger()`) and
any writes to `Base.stdout`/`Base.stderr` (by redirecting them to `devnull`). The macro
returns the value of `expr`.

Usage

  - `@silent some_expression` — runs `some_expression` with logging and standard
    output/stderr suppressed.

Returns

  - The result of evaluating `expr`.

Notes

  - Suppresses messages produced by logging macros such as `@info`, `@warn`, etc., by
    temporarily installing a `Logging.NullLogger()`.
  - Also redirects `Base.stdout` and `Base.stderr` to `devnull`, suppressing plain
    `print`/`println` output.
  - Evaluation happens in the caller's scope (macro hygiene preserved where appropriate);
    side effects of `expr` still occur, except for emitted text.
  - The suppression is lexically limited to the execution of `expr`; the original logger and
    standard streams are restored afterwards.
  - Exceptions raised by `expr` are propagated; suppression is still undone on exit.
"""
macro silent(expr)
    quote
        local _val
        # Mute logging macros (@info/@warn) no matter what global_logger was set to:
        with_logger(Logging.NullLogger()) do
            # Also mute anything writing to Base.stdout/Base.stderr:
            redirect_stdout(devnull) do
                redirect_stderr(devnull) do
                    _val = $(esc(expr))
                end
            end
        end
        _val
    end
end

export quiet_call, @silent, unmute_println, unmute_printf

end
