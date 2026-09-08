# ------------------------------------------------------------
# Metropolis-Hastings Performance 
# ------------------------------------------------------------
# This script evaluates the performance of the Metropolis-Hastings
# algorithm across different prior specifications and sample sizes.
# ------------------------------------------------------------

# install.packages("latex2exp")
library(latex2exp)


# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------

#' @description
#' Evaluates the performance of the Metropolis-Hastings algorithm for a given weight function.
#'
#' @details
#' Fits a Bayesian context tree model and runs a Metropolis-Hastings sampler using
#' `metropolis_vlmc()`. Computes the total posterior mass recovered by the chain,
#' the normalized error of the Maximum A Posteriori (MAP) tree probability, and the number
#' of unique sampled trees.
#'
#' @param w Prior weight function.
#' @param sample A sequence sample used for model fitting and sampling.
#'
#' @returns
#' A list containing:
#' - Total posterior probability mass recovered by the MH chain.
#' - Normalized difference (error) for the MAP tree probability.
#' - Number of unique trees visited by the MH chain.
#' 
mh_perf <- function(w, sample) {
  
  # Fit the Bayesian model
  model <- baConTree$new(sample, maximalDepth = 6, alpha = 0.5, priorWeights = w)
  
  # Obtain the MH chain
  set.seed(2)
  ch <- metropolis_vlmc(sample, num_samples = 10000, maximalDepth = 6, alpha = 0.5, priorWeights = w, burn_in = 1000)
  
  # Obtain the total posterior mass recovered by the MH chain
  prob_sum <- sum(sapply(ch$df$tree_contexts, function(contexts) {
    model$activateFromContexts(contexts)
    model$activeTreeProbabilities()[[2]]
  }))
  
  # Activate MAP tree
  model$activateMap()
  
  # Compute the normalized MAP error
  dif_map <- abs(model$activeTreeProbabilities()[[2]] - ch$df$prob[1]) / prob_sum
  
  return(list(
    prob_sum, 
    dif_map, 
    length(ch$df$tree_contexts)
  ))
}


#' @description
#' Plots benchmark metrics across sample sizes for four prior specifications.
#'
#' @param values_u Vector of metrics for the uniform prior U.
#' @param values_tcl Vector of metrics for the TCL prior.
#' @param values_ctw Vector of metrics for the CTW prior.
#' @param values_e Vector of metrics for the exponential prior.
#' @param ylab Character string for the y-axis label.
#' @param ylim Optional numeric vector of length 2 for y-axis limits.
plot_results <- function(values_u,
                         values_tcl,
                         values_ctw,
                         values_e,
                         ylab,
                         ylim = NULL) {
  
  plot(
    sample_sizes,
    values_u,
    type = "b",
    pch = 16,
    col = "blue",
    xaxt = "n",
    yaxt = ifelse(is.null(ylim), "s", "n"),
    xlab = TeX("$n$"),
    ylab = ylab,
    ylim = ylim
  )
  
  lines(sample_sizes, values_tcl, type = "b", pch = 17, col = "red")
  lines(sample_sizes, values_ctw, type = "b", pch = 15, col = "darkgreen")
  lines(sample_sizes, values_e,   type = "b", pch = 18, col = "orange")
  
  axis(1, at = sample_sizes, labels = sample_sizes)
  
  if (!is.null(ylim)) {
    axis(2, at = seq(ylim[1], ylim[2], length.out = 11))
  }
  
  legend(
    "bottomright",
    legend = c(TeX("$U$"), TeX("$T_2^3$"), TeX("$C$"), TeX("$E_{-1/2}$")),
    col = c("blue", "red", "darkgreen", "orange"),
    pch = c(16, 17, 15, 18),
    lty = 1,
    bty = "n"
  )
}


# ------------------------------------------------------------
# Simulation Setup & Data Generation
# ------------------------------------------------------------

# Sample sizes
sample_sizes <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)

# Generating samples
samples <- generate_samples(sample_sizes, true_contexts, true_probs)


# ------------------------------------------------------------
# Model Execution for Each Prior
# ------------------------------------------------------------

results_u   <- lapply(samples, function(sample) mh_perf(u, sample))
results_tcl <- lapply(samples, function(sample) mh_perf(tcl(3, 3), sample))
results_ctw <- lapply(samples, function(sample) mh_perf(ctw(6), sample))
results_e   <- lapply(samples, function(sample) mh_perf(e_alpha(0.5), sample))


# ------------------------------------------------------------
# Extract Metrics
# ------------------------------------------------------------

# Total probability
prob_u   <- sapply(results_u,   function(x) x[[1]])
prob_tcl <- sapply(results_tcl, function(x) x[[1]])
prob_ctw <- sapply(results_ctw, function(x) x[[1]])
prob_e   <- sapply(results_e,   function(x) x[[1]])

# Difference in MAP probabilities
dif_u   <- sapply(results_u,   function(x) x[[2]])
dif_tcl <- sapply(results_tcl, function(x) x[[2]])
dif_ctw <- sapply(results_ctw, function(x) x[[2]])
dif_e   <- sapply(results_e,   function(x) x[[2]])

# Number of sampled trees
n_u   <- sapply(results_u,   function(x) x[[3]])
n_tcl <- sapply(results_tcl, function(x) x[[3]])
n_ctw <- sapply(results_ctw, function(x) x[[3]])
n_e   <- sapply(results_e,   function(x) x[[3]])


# ------------------------------------------------------------
# Visualization & Results Summary
# ------------------------------------------------------------

# Plot total posterior probability mass
plot_results(
  prob_u,
  prob_tcl,
  prob_ctw,
  prob_e,
  ylab = TeX("$P_{mh}$"),
  ylim = c(0, 1)
)

# MAP probability error table
dif_map_table <- data.frame(
  sample_size = sample_sizes,
  u   = dif_u,
  td  = dif_tcl,
  ctw = dif_ctw,
  exp = dif_e
)

print(dif_map_table)

# Number of unique sampled trees table
n_trees_table <- data.frame(
  sample_size = sample_sizes,
  u   = n_u,
  td  = n_tcl,
  ctw = n_ctw,
  exp = n_e
)

print(n_trees_table)