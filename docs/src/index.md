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

`bld/`: Contains the generated files (after generating them). The folder is not part of the repository, but is created when running (certain parts of) the toolbox. That is, the folder contains:
- the generated files from the examples as subfolders of `bld/`.

`docs/`: Contains the documentation of the toolbox, that is:
- the source code in `src/`,
- and the generated documentation in `build/`.

`test/`: Contains the tests for the toolbox.

### Building the documentation

You can build the documentation *locally* by starting a new Julia REPL in the root directory of the repository, activating the environment, and running the following command: `include("docs/make.jl")`. You can access the documentation, once it is built locally, via running `python3 -m http.server --directory docs/build/`. If you then open your browser at [http://localhost:8000](http://localhost:8000), the documentation should render properly. Beyond that, the documentation is hosted via GitHub Pages and can be accessed [here](https://hildebrandecon.github.io/BASEtoolbox.jl/).

### Getting started with your model

The backbone of the toolbox is a computation algorithm to efficiently solve one- or two-asset heterogeneous agent models. The household problem, including all notation, is described in detail in [Household Problem](HouseholdProblem.md). The algorithm is described in [Computational Notes](ComputationalNotes.md).

If you want to add a new model, the recommended way is to start by copying one of the provided examples into a new folder in `examples/`. This way, you can make sure that all necessary files are present and that the toolbox can be run without any issues.

We provide a detailed description of the user inputs in a typical example in [General example structure](GeneralStructure.md). To decide which example is best suited as a starting point for your needs, you can look at the list of the provided examples in [Examples](GeneralStructure.md).

!!! tip
    If your model differs only in the aggregate model part, you can simply stay in the `examples/<your_example>/` folder. In this case you should not need to change the files `src/`. This also holds for some options built into the household problem already, e.g. taxes. Those can entirely be adjusted within your `example/<your_example>/` folder without changing the `src/` files.

!!! warning
    We only advise very proficient "HANK users" to change the `src/` files. This is because changing code at one place may not be enough for the model to be still valid.

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
