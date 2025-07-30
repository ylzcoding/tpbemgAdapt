#' call getsamples.full function
#' @export
fullGibbs = function(burnin_num, output_num, X, y,
                     a0, b0, beta0, nu0, xi0 = NULL, lambda0,
                     omega0, w0, sigmaSq0, omega_hyper,
                     tunea = 0.1, tuneb = 0.2,
                     uppera = 1, upperb = 1, 
                     woodbury = FALSE, approx = FALSE, diagX = FALSE, 
                     prshapeo, prrateo, prshapes, prrates){
  ### burn in and sampling
  samples = getsamples.full(num = output_num, X = X, y = y,
                            a0 = a0, b0 = b0,
                            beta0 = beta0, nu0 = nu0,
                            lambda0 = lambda0, omega0 = omega0,
                            w0 = w0, sigmaSq0 = sigmaSq0,
                            omega_hyper = omega_hyper, burn = burnin_num,
                            xi0 = xi0, tunea = tunea, tuneb = tuneb,
                            uppera = uppera, upperb = upperb, 
                            woodbury = woodbury, approx = approx, diagX = diagX, 
                            prshapeo = prshapeo, prrateo = prrateo,
                            prshapes = prshapes, prrates = prrates)
  return(samples)
}
