# Helper functions for TreeBH
# These functions are used internally by the different TreeBH implementations

#' Calculate Simes p-value
#' @param pvals Numeric vector of p-values
#' @return Simes combined p-value
#' @noRd
get_simes_p <- function(pvals) {
  # Remove NA values
  clean_pvals <- pvals[!is.na(pvals)]
  
  # If no valid p-values, return NA
  if (length(clean_pvals) == 0) {
    return(NA_real_)
  }
  
  # Sort p-values
  sorted_pvals <- sort(clean_pvals)
  n <- length(sorted_pvals)
  
  # Calculate Simes p-value: min_i( n * p_(i) / i )
  adjusted <- sorted_pvals * n / seq_along(sorted_pvals)
  min_val <- min(adjusted)
  
  return(min(1, min_val))
}

#' Calculate Fisher's combined p-value
#' @param pvals Numeric vector of p-values
#' @return Fisher's combined p-value
#' @noRd
get_fisher_p <- function(pvals) {
  # Remove NA values
  clean_pvals <- pvals[!is.na(pvals)]
  
  # If no valid p-values, return NA
  if (length(clean_pvals) == 0) {
    return(NA_real_)
  }
  
  # Calculate Fisher's statistic: -2 * sum(log(p))
  statistic <- -2 * sum(log(clean_pvals))
  
  # Return p-value from chi-squared distribution with 2k degrees of freedom
  return(pchisq(statistic, df = 2 * length(clean_pvals), lower.tail = FALSE))
}

