calculate_NSE <- function(observed, simulated) {
  mean_observed <- mean(observed)
  numerator <- sum((observed - simulated)^2)
  denominator <- sum((observed - mean_observed)^2)
  nse <- 1 - (numerator / denominator)
  return(nse)
}