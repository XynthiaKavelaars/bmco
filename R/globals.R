# Suppress R CMD check notes about "no visible binding"
# These variables are used in non-standard evaluation contexts
utils::globalVariables(c(
  "grp", "id", "x_var", "y1", "y2", 
  "Intercept", "grp_x", "cluster"
))