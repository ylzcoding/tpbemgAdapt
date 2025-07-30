#' get posterior samples with a, b sampled from Metropolis-Hasting (full Bayes)
#' @param burn numeric, number of burn-in samples
#' @param num numeric, number of samples to generate after burn-in for each hyperparameter
#' @param xi0 initial value for xi, the auxiliary variable, xi0 = NULL if not needed
#' @param tunea numeric from 0-1, step size of a in Metropolis-Hasting
#' @param tuneb numeric from 0-1, step size of b in Metropolis-Hasting
#' @param uppera upper bound for a, a ~ uniform(0, uppera)
#' @param upperb upper bound for b, b ~ uniform(0, upperb)
#' @return a list of posterior samples generated, including beta, nu, phi, w, sigmaSq, a, b, xi, as well as acceptance rates of a and b
#' @export
getsamples.full = function(num, X, y, a0, b0, beta0, nu0, lambda0, omega0, w0, sigmaSq0, omega_hyper,
                           burn = 0, xi0 = NULL,
                           tunea = 0.1, tuneb = 0.2, uppera = 1, upperb = 1, 
                           woodbury = FALSE, approx = FALSE, diagX = FALSE, 
                           prshapeo, prrateo, prshapes, prrates){
  aug <- !is.null(xi0)
  p = dim(X)[2]
  n = dim(X)[1]
  # num: number of samples to generate after burn-in
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
  omega_samples = rep(NA, num + burn)
  omega_samples[1] = omega0
  w_samples = rep(NA, num + burn)
  w_samples[1] = w0
  sigmaSq_samples = rep(NA, num + burn)
  sigmaSq_samples[1] = sigmaSq0
  a_samples = rep(NA, num + burn)
  a.acc = rep(NA, num + burn)
  a_samples[1] = a0
  b_samples = rep(NA, num + burn)
  b.acc = rep(NA, num + burn)
  b_samples[1] = b0

  for (i in 2:(burn + num)){
    beta_samples[i, ] = Gibbs_beta(X = X, y = y, omega = omega_samples[i-1], sigmaSq = sigmaSq_samples[i-1], nu = nu_samples[i-1, ], lambda = lambda_samples[i-1, ], 
                                   woodbury = woodbury, approx = approx, diagX = diagX)
    if (aug) {
      nu_samples[i, ] = aug_Gibbs_nu(p = p, omega = omega_samples[i-1], beta = beta_samples[i, ], lambda = lambda_samples[i-1, ], xi = xi_samples[i-1, ], a = a_samples[i-1])
      xi_samples[i, ] = Gibbs_xi(a = a_samples[i-1], p = p, nu = nu_samples[i, ])
    } else {
      nu_samples[i, ] = Naive_Gibbs_nu(a = a_samples[i-1], p = p, omega = omega_samples[i-1], beta = beta_samples[i, ], lambda = lambda_samples[i-1, ])
    }
    lambda_samples[i, ] = Gibbs_lambda(b = b_samples[i-1], p = p, omega = omega_samples[i-1], beta = beta_samples[i, ], nu = nu_samples[i, ])
    omega_samples[i] = Gibbs_phi(p, nu_samples[i, ], lambda_samples[i, ], w_samples[i-1], beta_samples[i, ], prshape = prshapeo, prrate = prrateo)
    w_samples[i] = Gibbs_w(omega_samples[i], omega_hyper)
    sigmaSq_samples[i] = Gibbs_sigmaSq(n, p, y, X, beta_samples[i, ], prshape = prshapes, prrate = prrates)
    if (tunea != 0) {
    a.mh = mh(p, a_samples[i-1], nu_samples[i, ], tunea, upper = uppera)
    a.acc[i] = a.mh$gen != a_samples[i-1]
    a_samples[i] = a.mh$gen
    } else {
      a.acc[i] = 0
      a_samples[i] = a_samples[i-1]
    }
    if (tuneb != 0) {
    b.mh = mh(p, b_samples[i-1], lambda_samples[i, ], tuneb, upper = upperb)
    b.acc[i] = b.mh$gen != b_samples[i-1]
    b_samples[i] = b.mh$gen
    } else {
      b.acc[i] = 0
      b_samples[i] = b_samples[i-1]
    }
  }
  res <- list(beta = beta_samples[(burn + 1):(burn + num), ],
              nu = nu_samples[(burn + 1):(burn + num), ],
              lambda = lambda_samples[(burn + 1):(burn + num), ],
              omega = omega_samples[(burn + 1):(burn + num)],
              w = w_samples[(burn + 1):(burn + num)],
              sigmaSq = sigmaSq_samples[(burn + 1):(burn + num)],
              a = a_samples[(burn + 1):(burn + num)],
              b = b_samples[(burn + 1):(burn + num)],
              a.acc = a.acc[(burn + 1):(burn + num)],
              b.acc = b.acc[(burn + 1):(burn + num)])
  if (aug) {
    res$xi = xi_samples[(burn + 1):(burn + num), ]
  }
  return(res)
}
