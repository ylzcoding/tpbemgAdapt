#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' This function automates the selection of starting points for the main EM_GibbsTPB
#' algorithm by running a competition between several candidate models (special cases
#' of the TPB prior). It selects the candidate that best fits the data and returns
#' its parameters as the starting points.
#'
#' @import glmnet
#' @param X, y Model matrix and response vector.
#' @param adapt_iter Max number of iterations for the mini-EM of the adaptive initialization
#' @param adapt_burnin Burn-in iterations for the "clean" Gibbs run for scoring.
#' @param adapt_samples Sampling iterations for the "clean" Gibbs run for scoring.
#' @param omega_init_guess A single numeric value for the initial guess of omega.
#' @param sigmaSq_init_guess A single numeric value for the initial guess of sigmaSq.
#' @param ... All other necessary parameters and functions to be passed to run_em_engine
#'   (e.g., IS_period, converge_rule, getsamples.emp, log_lik, M.step functions, etc.).
#' @return A list with the winning starting values: a.start, b.start, omega.start, sigmaSq.start
#' @export
initialize_adaptive <- function(X, y, 
                                omega_init_ebrr,
                                sigmaSq_init_ebrr,
                                iter_burnin, iter_samples,
                                adapt_burnin, adapt_samples,
                                adapt_iter, woodbury, IS_period = 1,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, converge_rule = "hyper",
                                approx = FALSE, option = "augmented", diagX = FALSE) {
  
  cat("--- Starting Data-Adaptive Initialization via Model Competition ---\n")
  n = nrow(X)
  
  # Collect all arguments to pass down to the engine
  # These are common to all engine calls within this initialization phase.
  engine_args <- list(
    X = X, y = y,
    max_iter = adapt_iter,
    iter_burnin = iter_burnin, # These are for the EM steps
    iter_samples = iter_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size, converge_rule = converge_rule,
    woodbury = woodbury, approx = approx,
    IS_period = IS_period, option = option, diagX = diagX
  )
  
  # List to store results from each candidate
  candidates <- list()
  
  # --- Candidate 1: Horseshoe (a=0.5, b=0.5) ---
  cat("Analyzing Candidate 1: Horseshoe...\n")
  hs_specific_args <- list(a_init = 0.5, b_init = 0.5,
                           omega_init = omega_init_ebrr, sigmaSq_init = sigmaSq_init_ebrr,
                           null.a = FALSE, null.b = FALSE, null.omega = TRUE, null.sigmaSq = TRUE)
  
  res_hs <- do.call(run_em_engine, c(hs_specific_args, engine_args))
  
  cat("Generating clean samples for Horseshoe scoring...\n")
  clean_samples_hs <- getsamples.emp(num = adapt_samples, X = X, y = y,
                                     a = res_hs$params$a, b = res_hs$params$b, 
                                     omega = res_hs$params$omega, sigmaSq = res_hs$params$sigmaSq,
                                     beta0 = res_hs$final_state$beta0, nu0 = res_hs$final_state$nu0, 
                                     lambda0 = res_hs$final_state$lambda0, burn = adapt_burnin,
                                     xi0 = res_hs$final_state$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
  
  beta_hs <- colMeans(clean_samples_hs$beta)
  score_hs <- get_initialization_score(X, y, beta_hs, res_hs$params$sigmaSq)
  
  candidates$horseshoe <- list(params = res_hs$params, score = score_hs)
  cat("Horseshoe Score:", score_hs, "\n")
  
  
  # --- Candidate 2: Student's t (approximated with large a, b=1) ---
  cat("Analyzing Candidate 2: Student's t (df=2)...\n")
  st_specific_args <- list(a_init = 100, b_init = 1.0, 
                           omega_init = omega_init_ebrr, sigmaSq_init = sigmaSq_init_ebrr,
                           null.a = FALSE, null.b = FALSE, null.omega = TRUE, null.sigmaSq = TRUE)
  res_st <- do.call(run_em_engine, c(st_specific_args, engine_args))
  
  cat("Generating clean samples for Student's t scoring...\n")
  clean_samples_st <- getsamples.emp(num = adapt_samples, X = X, y = y,
                                     a = res_st$params$a, b = res_st$params$b, 
                                     omega = res_st$params$omega, sigmaSq = res_st$params$sigmaSq,
                                     beta0 = res_st$final_state$beta0, nu0 = res_st$final_state$nu0, 
                                     lambda0 = res_st$final_state$lambda0, burn = adapt_burnin,
                                     xi0 = res_st$final_state$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
  
  beta_st <- colMeans(clean_samples_st$beta)
  score_st <- get_initialization_score(X, y, beta_st, res_st$params$sigmaSq)
  
  st_start_params <- res_st$params
  st_start_params$a <- 5.0 # Use a more moderate 'a' for the final start
  candidates$student_t <- list(params = st_start_params, score = score_st)
  cat("Student's t Score:", score_st, "\n")
  
  
  # --- Candidate 3: Normal-Gamma (approximated with large b, a=1) ---
  cat("Analyzing Candidate 3: Normal-Gamma...\n")
  ng_specific_args <- list(a_init = 1.0, b_init = 100, 
                           omega_init = omega_init_ebrr, sigmaSq_init = sigmaSq_init_ebrr,
                           null.a = FALSE, null.b = FALSE, null.omega = TRUE, null.sigmaSq = TRUE)
  res_ng <- do.call(run_em_engine, c(ng_specific_args, engine_args))
  
  cat("Generating clean samples for Normal-Gamma scoring...\n")
  clean_samples_ng <- getsamples.emp(num = adapt_samples, X = X, y = y,
                                     a = res_ng$params$a, b = res_ng$params$b, 
                                     omega = res_ng$params$omega, sigmaSq = res_ng$params$sigmaSq,
                                     beta0 = res_ng$final_state$beta0, nu0 = res_ng$final_state$nu0, 
                                     lambda0 = res_ng$final_state$lambda0, burn = adapt_burnin,
                                     xi0 = res_ng$final_state$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
  
  beta_ng <- colMeans(clean_samples_ng$beta)
  score_ng <- get_initialization_score(X, y, beta_ng, res_ng$params$sigmaSq)
  
  ng_start_params <- res_ng$params
  ng_start_params$b <- 5.0 # Use a more moderate 'b' for the final start
  candidates$normal_gamma <- list(params = ng_start_params, score = score_ng)
  cat("Normal-Gamma Score:", score_ng, "\n")
  
  
  # --- Candidate 4: Strawderman-Berger (a=1, b=0.5) ---
  cat("Analyzing Candidate 4: Strawderman-Berger...\n")
  sb_specific_args <- list(a_init = 1, b_init = 0.5,
                           omega_init = omega_init_ebrr, sigmaSq_init = sigmaSq_init_ebrr,
                           null.a = FALSE, null.b = FALSE, null.omega = TRUE, null.sigmaSq = TRUE)
  
  res_sb <- do.call(run_em_engine, c(sb_specific_args, engine_args))
  
  cat("Generating clean samples for Strawderman-Berger scoring...\n")
  clean_samples_sb <- getsamples.emp(num = adapt_samples, X = X, y = y,
                                     a = res_sb$params$a, b = res_sb$params$b, 
                                     omega = res_sb$params$omega, sigmaSq = res_sb$params$sigmaSq,
                                     beta0 = res_sb$final_state$beta0, nu0 = res_sb$final_state$nu0, 
                                     lambda0 = res_sb$final_state$lambda0, burn = adapt_burnin,
                                     xi0 = res_sb$final_state$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
  
  beta_sb <- colMeans(clean_samples_sb$beta)
  score_sb <- get_initialization_score(X, y, beta_sb, res_sb$params$sigmaSq)
  
  candidates$strawderman_berger <- list(params = res_sb$params, score = score_sb)
  cat("Strawderman-Berger Score:", score_sb, "\n")
  
  
  # --- Candidate 5: Inv_Gamma_Gamma (a=0.5+1/n, b=1/n) ---
  cat("Analyzing Candidate 5: Inv_Gamma_Gamma...\n")
  igg_specific_args <- list(a_init = 0.5+1/n, b_init = 1/n,
                           omega_init = omega_init_ebrr, sigmaSq_init = sigmaSq_init_ebrr,
                           null.a = FALSE, null.b = FALSE, null.omega = TRUE, null.sigmaSq = TRUE)
  
  res_igg <- do.call(run_em_engine, c(igg_specific_args, engine_args))
  
  cat("Generating clean samples for Inv_Gamma_Gamma scoring...\n")
  clean_samples_igg <- getsamples.emp(num = adapt_samples, X = X, y = y,
                                     a = res_igg$params$a, b = res_igg$params$b, 
                                     omega = res_igg$params$omega, sigmaSq = res_igg$params$sigmaSq,
                                     beta0 = res_igg$final_state$beta0, nu0 = res_igg$final_state$nu0, 
                                     lambda0 = res_igg$final_state$lambda0, burn = adapt_burnin,
                                     xi0 = res_igg$final_state$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
  
  beta_igg <- colMeans(clean_samples_igg$beta)
  score_igg <- get_initialization_score(X, y, beta_igg, res_igg$params$sigmaSq)
  
  candidates$inv_gamma_gamma <- list(params = res_igg$params, score = score_igg)
  cat("Inv_Gamma_Gamma Score:", score_igg, "\n")
  
  
  
  ## --- Candidate 6: Ridge Regression (a, b -> inf)---
  #cat("Analyzing Candidate 4: Ridge Regression...\n")
  
  #cv_ridge <- cv.glmnet(X, y, alpha = 0, intercept = FALSE)
  #best_lambda <- cv_ridge$lambda.min
  #final_ridge_model <- glmnet(X, y, alpha = 0, lambda = best_lambda, intercept = FALSE)
  #beta_ridge <- as.numeric(coef(final_ridge_model)[-1])
  #y_hat_ridge <- predict(final_ridge_model, s = best_lambda, newx = X)
  #sigmaSq_ridge <- max(1e-6, var(y - y_hat_ridge))
    
  ## Heuristic for omega
  #omega_ridge <- mean(beta_ridge^2) / sigmaSq_ridge 
  #if (omega_ridge < 1e-8) omega_ridge <- 1e-8 # Prevent numerical issues
  #score_ridge <- get_initialization_score(X, y, beta_ridge, sigmaSq_ridge)
    
  ## For the final start, use more moderate values for a and b
  #ridge_start_params <- list(a.start = 5.0, b.start = 5.0, 
                             #omega.start = omega_ridge, sigmaSq.start = sigmaSq_ridge)
    
  #candidates$ridge <- list(params = ridge_start_params, score = score_ridge)
  #cat("Ridge Score:", score_ridge, "\n")
  
  
  # --- Competition and Decision ---
  
  scores <- sapply(candidates, `[[`, "score")
  winner_name <- names(which.max(scores))
  winning_params <- candidates[[winner_name]]$params
  cat("--- Initialization Complete. Winner is:", winner_name, "---\n")
  
  return(
    list(winning_params = winning_params,
         winner_name = winner_name)
  )
}