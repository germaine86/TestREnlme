## Suppress R CMD check NOTEs for ggplot2 non-standard evaluation
## column names used inside aes(); these are data.frame columns
## constructed locally within each plotting function, not undefined
## globals.
utils::globalVariables(c(
  "subject", "kappa_val", "status", "value", "type",
  "Observed", "Fitted", "fitted", "residual", "parameter",
  "estimate", "method", "m"
))
