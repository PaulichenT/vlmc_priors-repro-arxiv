# ------------------------------------------------------------
# Simulation study — Scenario (b)
# ------------------------------------------------------------
# This script contains all simulation results concerning
# scenario (b).
# ------------------------------------------------------------


# ------------------------------------------------------------
# Generator model (VLMC model (b))
# ------------------------------------------------------------

# True contexts
true_contexts <- c(
  "*.0",
  "*.1.0",
  "*.1.1.0",
  "*.1.1.1.0",
  "*.1.1.1.1"
)

# Transition probabilities
true_probs <- list(
  c(0.1, 0.9),
  c(0.5, 0.5),
  c(0.5, 0.5),
  c(0.5, 0.5),
  c(0.9, 0.1)
)


# ------------------------------------------------------------
# Prior specifications (context-tree functions `f`)
# ------------------------------------------------------------

priors <- list(
  LDEP4  = d_l(4),
  CTW_05 = ctw(10),
  BCT_02 = b_beta(0.2, 10, 2),
  BCT_07 = b_beta(0.7, 10, 2),
  TDD_3  = tcl(3, 2),
  I_0    = i_a("0"),
  TDD_4  = tcl(4, 2),
  EXP_2  = e_alpha(2),
  EXP_5  = e_alpha(5),
  EXP_N  = e_ls
)


# ------------------------------------------------------------
# Sample sizes
# ------------------------------------------------------------

sample_sizes <- c(200, 500, 1000, 2500)


# ------------------------------------------------------------
# Generate samples
# ------------------------------------------------------------

samples <- generate_samples(sample_sizes, true_contexts, true_probs)


# ------------------------------------------------------------
# Evaluation metrics (Table 2)
# ------------------------------------------------------------

results_df <- run_scenario_from_samples(samples, sample_sizes, priors, true_contexts)
print(results_df)


# ------------------------------------------------------------
# Maximal depth selection
# ------------------------------------------------------------

selected_depths <- lapply(samples, select_depth_bf)
print(selected_depths)


# ------------------------------------------------------------
# Model selection
# ------------------------------------------------------------

best_models <- lapply(samples, function(sample) {model_selection(sample, priors = priors)})
print(best_models)

