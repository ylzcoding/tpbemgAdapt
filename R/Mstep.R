#' E-M steps
#' @import bzinb
#' @param sample_matrix a matrix with dimension num_sample*p - nu for estimating a and lambda for estimating b
#' @export
M.step = function(sample_matrix){
  # sample_matrix: a num_sample*p matrix
  empMean_log = colMeans(log(sample_matrix))
  return(bzinb::idigamma(mean(empMean_log)))
}

#' @param beta_matrix beta matrix with dimension num_sample*p
#' @export
M.step_sigmaSq = function(beta_matrix, X, y, n, diagX = FALSE){
  
  beta_samples = t(beta_matrix)
  num_samples = dim(beta_samples)[2]
  if (!diagX) {
    beta_beta_prime = vector("list", length = num_samples)
    for(i in 1:num_samples){
      beta_beta_prime[[i]] = beta_samples[, i] %*% t(beta_samples[, i])
    }
    empMean_beta_beta_prime = Reduce("+", beta_beta_prime)/ num_samples
    empMean_beta = colMeans(beta_matrix)
    res = (sum(diag(t(X) %*% X %*% empMean_beta_beta_prime)) - 2 * t(y) %*% X %*% empMean_beta + t(y) %*% y)/ n
  } else {
    empMean_beta = colMeans(beta_matrix)
    empMean_beta_beta_prime = colMeans(beta_matrix^2)
    res = (sum(diag(X)^2*empMean_beta_beta_prime) - 2 * sum(y*diag(X)*empMean_beta) + sum(y^2))/ n
  }
  return(max(res, 1e-16))
}

#' @param sigmaSq numeric, estimated sigmaSq from last iteration
#' @export
M.step_phi = function(beta_matrix, lambda_matrix, nu_matrix){
  empMean = colMeans(beta_matrix^2 * lambda_matrix/ nu_matrix)
  return(mean(empMean))
}
