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
# Prior specifications (weight functions `w`)
# ------------------------------------------------------------

priors <- list(
  LDEP4  = d_l_m(0, 4),
  CTW    = ctw(10),
  BCT_02 = bct_beta(0.2, 10, 2),
  BCT_07 = bct_beta(0.7, 10, 2),
  TDD_3  = t_beta_l(3, 2),
  I_0    = i_a("0"),
  TDD_4  = t_beta_l(4, 2),
  EXP_2  = k_beta(exp(-2)),
  EXP_5  = k_beta(exp(-5)),
  EXP_N  = e_beta(-1)
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
# Maximal depth selection (Section A3, Supplementary material)
# ------------------------------------------------------------

selected_depths <- lapply(samples, select_depth)
print(selected_depths)


