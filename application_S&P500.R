# ------------------------------------------------------------
# Application - S&P500
# ------------------------------------------------------------
# This script contains the results of the application to
# financial markets.
# ------------------------------------------------------------


# ------------------------------------------------------------
# Data sequences
# ------------------------------------------------------------

sequence_return <- readRDS("sequence_return.rds")
sequence_volatility <- readRDS("sequence_volatility.rds")


# ------------------------------------------------------------
# Prior specifications (weight functions `w`)
# ------------------------------------------------------------

priors <- list(
  CTW   = ctw(10),
  BCT25 = bct_beta(0.25, 10, 3), 
  BCT75 = bct_beta(0.75, 10, 3), 
  BCT90 = bct_beta(0.90, 10, 3),
  U     = u, 
  EXP   = e_beta(-0.2),
  TD33  = t_beta_l(3, 3)
)


# ------------------------------------------------------------
# Results for the daily return data
# ------------------------------------------------------------

ret_values <- log_evidence(sequence_return, 10, 0.5, priors)
print(ret_values)


# ------------------------------------------------------------
# Results for the daily volatility data
# ------------------------------------------------------------

vol_values <- log_evidence(sequence_volatility, 10, 0.5, priors)
print(vol_values)
