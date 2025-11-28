adaTPB: Adaptive Normal Three-Parameter Beta Prior for High-Dimensional Regression

adaTPB implements a robust Empirical Bayes framework for high-dimensional linear regression ($p \gg n$). It utilizes the normal Three-Parameter Beta (TPB) prior with a novel model-based warm-start strategy to adaptively estimate hyperparameters.


Quick Start

Here is a minimal example of how to use the main function:
```
library(tpbemg)
```
1. Generate synthetic data

```
n <- 50 
p <- 200 
X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
beta_true <- c(rep(2, 5), rep(0, p - 5)) 
y <- X %*% beta_true + rnorm(n)
```


2. Fit the model using the adaptive TPB prior

The function automatically handles initialization and hyperparameter selection
```
fit <- EM_GibbsTPB_adapt(X = X, y = y,
                         iter_pre_opt = 500, pre_opt_burnin = 1000, pre_opt_samples = 1000,
                         iter_burnin_selection = 0, iter_selection = 5000,
                         mainEM_burnin = 1000, mainEM_samples = num_tpb_em,
                         output_burnin = num_burnin, output_samples = num_samples,
                         mainEM_iter = 500, 
                         IS_period = 1, option = "augmented",
                         woodbury = woodbury, approx = FALSE,
                         delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-4,
                         window_size = 5, converge_rule = "hyper",
                         a = NULL, b = NULL, omega = NULL, sigmaSq = NULL,
                         a.start = 0.5, b.start = 0.5, omega.start = omega_init_guess, sigmaSq.start = sigmaSq_init_guess,
                         adapt_init = TRUE, diagX = FALSE)
```

3. Inspect results
```
beta_hat <- colMeans(fit$beta)
plot(beta_hat, main = "Estimated Coefficients")

a_hat <- fit$final_params$a
b_hat <- fit$final_params$b
omega_hat <- fit$final_params$omega
```

Key Features

Adaptive: Automatically selects hyperparameters ($a, b, \omega, \sigma^2$) based on data structure.

Scalable: Efficient computation for high-dimensional data ($p \gg n$).

Robust: Features a model-based warm-start strategy to ensure convergence.
