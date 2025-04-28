# BASEforHANK.jl Documentation

## Introduction

This manual documents the Julia package **BASEforHANK**, which provides a toolbox for the Bayesian Solution and Estimation (BASE) of a heterogeneous-agent New-Keynesian (HANK) model.

It comes with examples that showcase how to use the package. Originally, the code accompanied the paper [Bayer, Born, and Luetticke (2024, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20201875). Note that the toolbox is not a 1-for-1 replication package for the linked paper. In particular, the preset resolution is smaller.

## First steps

### Installation

We recommend to use [Julia for VSCode IDE](https://www.julia-vscode.org) as a front-end to Julia. To get started with the toolbox, simply download or clone the folder, e.g. via `git clone` and set your `cd` to the project directory. Then start the Julia REPL and type `]` so that you can call
```julia-repl
(v1.10) pkg> activate .

(BASEtoolbox) pkg> instantiate
```
This will install all needed packages. For more on Julia environments, see [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project).

!!! warning
    Before you activate the environment, make sure that you are in the main directory, in which the `Project.toml` file is located. In case you accidentally activated the environment in a subfolder, empty `.toml` files will be created that you need to delete before proceeding in the correct folder.

We have tested the module on Julia Version 1.10.3 (macOS, `arm64-apple-darwin22.4.0`), Version 1.10.4 (Windows, `x86_64-w64-mingw32`), and Version 1.10.5 (Linux, `x86_64-linux-gnu`). You can find out which version you are using by calling `versioninfo()` in the Julia REPL.

If you use a different editor, make sure that the environment is correctly set, as otherwise the instantiated packages might not be found.

### Toolbox folder structure

In the following, we call the root directory of the repository `BASEtoolbox.jl` (which is the directory containing, for instance, `Project.toml`). The folder structure is as follows:

`src/`: Contains the source code of the toolbox, that is:
- the main module file `BASEforHANK.jl`,
- the submodules in the folder `SubModules/`,
- and the pre-processing functions in the folder `Preprocessor/`.

`examples/`: Contains the examples that showcase the toolbox. For each example, there is a subfolder in `examples/` that contains the main file to run the example as well as all relevant files for the example. The baseline example that showcases most functions of the toolbox is given by `examples/baseline/main.jl`. This is strictly required as it serves as the baseline for testing and documentation.

`bld`: Contains the generated files (after generating them). The folder is not part of the repository, but is created when running (certain parts of) the toolbox. That is, the folder contains:
- the generated files from the examples as subfolders of `bld/`.

`docs/`: Contains the documentation of the toolbox, that is:
- the source code in `src/`,
- and the generated documentation in `build/`.

### Building the documentation

You can build the documentation *locally* by starting a new Julia REPL in the root directory of the repository, activating the environment, and running the following command: `include("docs/make.jl")`. You can access the documentation, once it is built locally, via running `python3 -m http.server --directory docs/build/`. If you then open your browser at [http://localhost:8000](http://localhost:8000), the documentation should render properly. Beyond that, the documentation is hosted via GitHub Pages and can be accessed [here](https://hildebrandecon.github.io/BASEtoolbox.jl/).

## Setting up your model

### Getting started

If you want to add a new model, the recommended way is to start by copying one of the provided examples into a new folder in `examples/`. This way, you can make sure that all necessary files are present and that the toolbox can be run without any issues.

[to do: add more details here]

### Special features

#### Multiple economies or sectors

The structure of the problem at hand allows for multiple economies or sectors in the aggregate model part as it would be nothing else than stating additional equations for each of them. To simplify writing up such models, the preprocessor allows using a certain syntax to mark the equations that would need to be repeated or copied for each economy and numbered accordingly. Instead, the preprocessor _automatically_ copies marked equations. Because many ASCII characters are reserved for other purposes, the syntax uses the `?` symbol.

The HANK economy is considered to be economy 1, not receiving an additional number. Any repeated equation will then get an additional number, starting with 2. The preprocessor will then automatically replace the `?` symbol with the respective number and will repeat the equation accordingly.

E.g., the following two code snippets from `input_aggregate_model.mod` are equivalent for the preprocessor:

```julia
# TFP HANK economy
F[indexes.Z]    = (log(ZPrime)) - (m_par.ρ_Z * log(Z))

# TFP economy 2
F[indexes.Z2]   = (log(Z2Prime)) - (m_par.ρ_Z * log(Z2))
```
```julia
# TFP HANK and all other economies
F[indexes.Z?]    = (log(Z?Prime)) - (m_par.ρ_Z * log(Z?))
```

At the beginning of the `input_aggregate_model.mod` file, the user has to specify the number of economies or sectors in their model by the following line:
```julia
# Setting the number of economies to 2.
@R?2
```

The same logic applies when specifying the variable names in `input_aggregate_names.jl` and `input_aggregate_steady_state.mod`. Consequently, the syntax can be also used when referring to steady state variables during writing down the model. E.g., ```Z?SS``` refers to the steady state value of variable `Z` for each economy.

In order to enable convenient syntax highlighting and because the `?` character serves other purposes in `Julia`, `.mod` files are used instead of `.jl` files for the user-specified model equations. In VS Code, you can still set the syntax highlighting to Julia for `.mod` files by clicking on the language mode in the bottom right corner of the editor.

## Methods

The following gives an overview on the methods provided by the toolbox. They naturally build on each other. For a more detailed documentation of each method, including its main functions, see the respective section.

### Preprocessing, parsing, and incomes

[to do: add more details here]

### Steady state and preparing linearization

[to do: add more details here]

Given model and parameters, this method computes the steady state of the model, i.e. the stationary equilibrium. It is summarized in a `SteadyStateStruct`.

Based on the steady state, the linearization of the model is prepared, which involves dimensionality reduction. The output of this step is collected in a `struct` `SteadyResults`.

For more details, see [Steady State](SteadyState.md).

### Linearization and model reduction

[to do: add more details here]

### Estimation

[to do: add more details here]

### Plots and statistics

[to do: add more details here]
