# ------------------------------------------------------------
# Simulation study — Scenario (a)
# ------------------------------------------------------------
# This script contains all simulation results concerning
# scenario (a).
# ------------------------------------------------------------


# ------------------------------------------------------------
# Generator model (VLMC model (a))
# ------------------------------------------------------------

# True contexts
true_contexts <- c(
  "*.0.0.0",
  "*.0.0.1",
  "*.0.1.0",
  "*.0.1.1",
  "*.1.0.0",
  "*.1.0.1",
  "*.1.1"
)

# Transition probabilities
true_probs <- list(
  c(0.9, 0.1),
  c(0.6, 0.4),
  c(0.7, 0.3),
  c(0.3, 0.7),
  c(0.8, 0.2),
  c(0.4, 0.6),
  c(0.4, 0.6)
)


# ------------------------------------------------------------
# Prior specifications (context-tree functions `f`)
# ------------------------------------------------------------

priors <- list(
  LDEP3  = d_l(3),
  CTW    = ctw(10),
  BCT_02 = b_beta(0.2, 10, 2),
  BCT_07 = b_beta(0.7, 10, 2),
  TDD_32 = tcl(3, 2),
  TDD_38 = tcl(3, 8),
  TDD_33 = tcl(3, 3),
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
# Evaluation metrics (Table 3)
# ------------------------------------------------------------

results_df <- run_scenario_from_samples(samples, sample_sizes, priors, true_contexts)
print(results_df)


# ------------------------------------------------------------
# Maximal Depth selection
# ------------------------------------------------------------

selected_depths <- lapply(samples, select_depth_bf)
print(selected_depths)


# ------------------------------------------------------------
# Model selection
# ------------------------------------------------------------

best_models <- lapply(samples, function(sample) {model_selection(sample, priors = priors)})
print(best_models)
