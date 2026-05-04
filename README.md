# vlmc_priors-repro-arxiv

Simulation results associated with the paper [arXiv:2603.25806](https://arxiv.org/abs/2603.25806)

---

## Instructions

To reproduce the simulation results presented in the paper, follow the steps below.

### 1. Install the `bacontrees` package

Install the development version of the package from [GitHub](https://github.com/Freguglia/bacontrees) with:

``` r
# install.packages("pak")
pak::pak("Freguglia/bacontrees")
```

### 2. Compile auxiliary scripts

Run the following scripts **in order**, as they define required functions and dependencies:

1. `context-tree_functions.R`
2. `simulation_functions.R`
3. `model_selection_bf.R`

### 3. Run simulations

Execute the script corresponding to the desired simulation scenario. 

For example, to obtain the results for scenario (a), run `simulation_scenario_a.R`. 

