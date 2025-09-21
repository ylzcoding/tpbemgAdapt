#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' This function automates the selection of starting points for the main EM_GibbsTPB
#' algorithm by running a competition between several candidate models (special cases
#' of the TPB prior). It selects the candidate that best fits the data and returns
#' its parameters as the starting points.
#' 
#' @param X, y Model matrix and response vector.
#' @param omega_init_guess, sigmaSq_init_guess Initial guesses for omega and sigmaSq, typically derived from a fast Ridge fit.
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param pre_opt_burnin, pre_opt_samples burn-in and sample size for the Gibbs sampler within each pre-optimization EM step.
#' @param iter_burnin_selection Number of burn-in iterations for the main model selection MCMC chain.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param candidates A list defining the candidate models and their fixed (a, b) values.
#' @param ... Other arguments to be passed to run_em_engine and getsamples.emp.
#' @return A list containing the winning parameters and the name of the winning model.
#' @export

initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt = 200,
                                pre_opt_burnin = 1000, pre_opt_samples = 1000,
                                iter_burnin_selection = 0, # if we need some burn-in iterations for the model selection procedure
                                iter_selection = 1000,
                                candidates = list(
                                  horseshoe = list(a = 0.5, b = 0.5),
                                  sb = list(a = 0.5, b = 0.1),
                                  student_t = list(a = 10.0, b = 1.0),
                                  normal_gamma = list(a = 1.0, b = 10.0)
                                ),
                                woodbury = TRUE, IS_period = 1,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, converge_rule = "hyper",
                                approx = FALSE, option = "augmented", diagX = FALSE) {
  
  
  cat("Stage 1: Pre-optimizing each candidate model...\n")
  
  pre_optimized_params <- list()
  pre_optimized_states <- list()
  
  engine_args <- list(
    X = X, y = y, 
    max_iter = iter_pre_opt, 
    iter_burnin = pre_opt_burnin,
    iter_samples = pre_opt_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size, converge_rule = converge_rule,
    woodbury = woodbury, approx = approx,
    IS_period = IS_period, option = option, diagX = diagX)
  
  for (name in names(candidates)) {
    cat("Optimizing candidate model:", name, "...\n")
    candidate_ab <- candidates[[name]]
    specific_args <- list(
      a_init = candidate_ab$a, b_init = candidate_ab$b,
      omega_init = omega_init_guess, sigmaSq_init = sigmaSq_init_guess,
      null.a = FALSE, null.b = FALSE,
      null.omega = TRUE, null.sigmaSq = TRUE
    )
    opt_res <- do.call(run_em_engine, c(specific_args, engine_args))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
  }
  
  cat("Stage 2: Starting iterative model selection...\n")
  
  num_candidates <- length(candidates)
  model_names <- names(candidates)
  model_choices <- integer(iter_selection) # Store only post-burn-in choices
  
  current_model_idx <- sample(1:num_candidates, 1)
  
  # Get the first single sample
  current_states <- pre_optimized_states[[current_model_idx]]
  current_beta <- current_states$beta0
  current_nu <- current_states$nu0
  current_lambda <- current_states$lambda0
  current_xi <- current_states$xi0
  
  iter_total <- iter_burnin_selection + iter_selection
  
  for (j in 1:iter_total) {
    # Update model indicator variable, delta
    log_weights <- sapply(model_names, function(name) {
      params <- pre_optimized_params[[name]]
      calculate_complete_loglik(
        beta_vec = current_beta, nu_vec = current_nu, lambda_vec = current_lambda,
        X = X, y = y, model_params = params
      )
    })
    
    max_log_weight <- max(log_weights)
    weights <- exp(log_weights - max_log_weight)
    probs <- weights / sum(weights)
    
    current_model_idx <- sample(1:num_candidates, 1, prob = probs)
    
    # Update beta, nu, and lambda based on the newly selected model
    current_params <- pre_optimized_params[[current_model_idx]]
    new_sample <-  getsamples.emp(num = 1, X = X, y = y,
                                  a = current_params$a, b = current_params$b, 
                                  omega = current_params$omega, sigmaSq = current_params$sigmaSq,
                                  beta0 = current_beta, nu0 = current_nu, 
                                  lambda0 = current_lambda, burn = 1,
                                  xi0 = current_xi, woodbury=woodbury, approx=approx, diagX=diagX)
    
    current_beta <- new_sample$beta
    current_nu <- new_sample$nu
    current_lambda <- new_sample$lambda
    current_xi <- new_sample$xi
    
    # Only store choices after the burn-in period
    if (j > iter_burnin_selection) {
      storage_idx <- j - iter_burnin_selection
      model_choices[storage_idx] <- current_model_idx
    }
  }
  
  model_counts <- table(factor(model_choices, levels = 1:num_candidates, labels = model_names))
  winner_name <- names(which.max(model_counts))
  winning_params <- pre_optimized_params[[winner_name]]
  
  return(
    list(winning_params = winning_params,
         winner_name = winner_name)
  )
}
  
  
  
  
  
  
  
                                  