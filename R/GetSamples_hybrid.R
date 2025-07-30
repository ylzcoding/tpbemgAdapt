#' @param burn numeric, number of burn-in samples
#' @param num numeric, number of samples to generate after burn-in for each hyperparameter
#' @param xi0 initial value for xi, the auxiliary variable, xi0 = NULL if not needed
#' @return a list of posterior samples generated given a and b
#' @export
getsamples.hybrid = function(num, X, y, a, b, phi_hyper, phi0, w0, sigmaSq0, beta0, nu0, lambda0, burn, xi0 = NULL, diagX = FALSE){
  aug <- !is.null(xi0)
  p = dim(X)[2]
  n = dim(X)[1]

  beta_samples = matrix(NA, nrow = num + burn, ncol = p)
  beta_samples[1, ] = beta0
  nu_samples = matrix(NA, nrow = num + burn, ncol = p)
  nu_samples[1, ] = nu0
  if (aug) {
    xi_samples = matrix(NA, nrow = num + burn, ncol = p)
    xi_samples[1, ] = xi0
  }
  lambda_samples = matrix(NA, nrow = num + burn, ncol = p)
  lambda_samples[1, ] = lambda0
  phi_samples = rep(NA, num + burn)
  phi_samples[1] = phi0
  w_samples = rep(NA, num + burn)
  w_samples[1] = w0
  sigmaSq_samples = rep(NA, num + burn)
  sigmaSq_samples[1] = sigmaSq0

  for (i in 2:(num + burn)){
    beta_samples[i, ] = Gibbs_beta(X, y, phi_samples[i-1], sigmaSq_samples[i-1], nu_samples[i-1, ], lambda_samples[i-1, ],
                                  diagX = diagX)
    if (aug) {
      nu_samples[i, ] = aug_Gibbs_nu(p, phi_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], lambda_samples[i-1, ], xi_samples[i-1, ])
      xi_samples[i, ] = Gibbs_xi(a, p, nu_samples[i, ])
    } else {
      nu_samples[i, ] = Naive_Gibbs_nu(a, p, phi_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], lambda_samples[i-1, ])
    }
    lambda_samples[i, ] = Gibbs_lambda(b, p, phi_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], nu_samples[i, ])
    phi_samples[i] = Gibbs_phi(p, nu_samples[i, ], lambda_samples[i, ], w_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ])
    w_samples[i] = Gibbs_w(phi_samples[i], phi_hyper)
    sigmaSq_samples[i] = Gibbs_sigmaSq(n, p, y, X, beta_samples[i, ], phi_samples[i], nu_samples[i, ], lambda_samples[i, ])
  }
  res <- list(beta = beta_samples[(burn + 1):(burn + num), ],
              nu = nu_samples[(burn + 1):(burn + num), ],
              lambda = lambda_samples[(burn + 1):(burn + num), ],
              phi = phi_samples[(burn + 1):(burn + num)],
              w = w_samples[(burn + 1):(burn + num)],
              sigmaSq = sigmaSq_samples[(burn + 1):(burn + num)])
  if (aug) {
    res$xi = xi_samples[(burn + 1):(burn + num), ]
  }
  return(res)
}




