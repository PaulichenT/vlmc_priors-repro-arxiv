# `vlmc_priors-repro-arxiv`

Computational code and results associated with the preprint [*Context Tree Prior Distributions Based on Node Weighting with Exact Bayes Factors*](https://arxiv.org/abs/2603.25806).

---

## Reproducing the Results

Follow the steps below to reproduce the simulation results, real-data application, and supplementary analyses presented in the paper.

### 1. Install the `bacontrees` package

Install the development version of the [`bacontrees`](https://github.com/Freguglia/bacontrees) package from GitHub:

```r
# Install pak, if necessary
install.packages("pak")

# Install the development version of bacontrees
pak::pak("Freguglia/bacontrees")
```

### 2. Compile the auxiliary scripts

Before running the analyses, execute the following scripts **in the specified order**:

1. `weight_functions.R`
2. `simulation_functions.R`
3. `model_selection.R`

These scripts define the weight functions, simulation utilities, and model-selection procedures required by the subsequent analyses.

### 3. Main analyses and results

#### Simulation study

To reproduce the results for a particular simulation scenario, run the corresponding script.

For example, to reproduce **Scenario (a)**:

`simulation_scenario_a.R`

Similarly, run the script corresponding to the desired scenario.

#### Real-data application

To reproduce the analysis of the **S&P 500** data, run:

`application_S&P500.R`

This script executes the complete application pipeline and produces the results reported in the paper.

### 4. Supplementary analysis

To evaluate the performance of the **Metropolis–Hastings algorithm**, run:

`MH_perf.R`.



