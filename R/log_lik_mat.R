#' @param beta_samples S*p matrix, posterior beta samples
#' @param sigmaSq numeric, estimated sigmaSq
#' @return log-likelihood matrix, (i, j)-entry is the log-likelihood of y_i under the j-th sampled beta.
#' @export
log_lik <- function(X, y, beta_samples, sigmaSq) {
  pred_mat <- X %*% t(beta_samples)
  y_mat <- matrix(rep(y, ncol(pred_mat)), nrow = nrow(X)) # y_mat: n*S matrix, each column is a replication of y
  log_lik_mat <- dnorm(y_mat, mean = pred_mat, sd = sqrt(sigmaSq), log = TRUE)
  return(log_lik_mat)
}

#' Wrapper function to calculate the score for a single initialization.
#' @param beta_k posterior mean after a given initialization and moderate number of iterations.
#' @param sigmaSq_k estimated sigmaSq after a given initialization and moderate number of iterations.
#' @return numeric, observed-data log-likelihood
#' @export
get_initialization_score <- function(X, y, beta_k, sigmaSq_k) {
  # Converts a point estimate beta_k into a 1-row matrix to be used
  # by the user's log_lik function, then sums the result to get a single score.
  
  beta_k_mtx <- matrix(beta_k, nrow = 1)
  log_lik_matrix <- log_lik(X = X, y = y, beta_samples = beta_k_mtx, sigmaSq = sigmaSq_k)
  log_lik_score <- sum(log_lik_matrix)
  return(log_lik_score)
}
