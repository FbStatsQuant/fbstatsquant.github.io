library(ggplot2)
library(tidyr)
library(dplyr)
library(Matrix)
library(Rcpp)
library(zoo)

# source("Sequential_Joint_Quantiles_Scale.R")


cat("Compiling C++ ADMM solver...\n")
Rcpp::sourceCpp("C:/Users/felip/OneDrive - Rice University/Proximal MCMC/Project/Code/Core/prox_admm_optimized_fixed.cpp")
cat("C++ compilation successful!\n\n")

set.seed(100)
n_data <- 400
x <- seq(1, n_data, length.out = n_data)
x_idx  <- x
xs <- x / n_data

# ---------- Model / noise selection (mirrors Individual_Sampler_Scale.R) ----------
model      <- 1          # 1 or 2
noise_name <- "Beta"   # "Normal" / "Beta" / "Mixed"

# Build baseline + lambda prior parameters per model (same logic as individual)
if (model == 1) {
  m1 <- n_data / 100
  baseline <- 2.5 * (x_idx >= 1        & x_idx <= m1*20)  +
              1.0 * (x_idx >= m1*20+1  & x_idx <= m1*40)  +
              3.5 * (x_idx >= m1*40+1  & x_idx <= m1*60)  +
              1.5 * (x_idx >= m1*60+1  & x_idx <= m1*100)
  C_1   <- 1.0
  #a_lam <- C_1 * n_data
  #b_lam <- n_data
  order <- 0
  a_lam <- 1
  b_lam <- C_1*n_data^(2/3)
  
} else if (model == 2) {
  C_1      <- 1
  baseline <- 2 + sin(4 * xs - 2) + 2 * exp(-30 * (4 * xs - 2)^2)
  #a_lam <- 1
  #b_lam <- (1/c_1) * sqrt(n_data)
  order <- 1
  a_lam <- C_1*n_data
  b_lam <- n_data^(1/2)
} else if (model == 3){
  C_1 <- 1
  baseline <- sin(2 * pi * xs)
  a_lam <- C_1
  b_lam <- 1
  order <- 2
  
}


# ---------- Quantile grid ----------
tau_vec <- c(0.002, 0.004, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
K <- length(tau_vec)

# Upper quant tau_vec <- c(0.3, 0.5, 0.7, 0.8, 0.9, 0.925, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.993, 0.996, 0.999)

# Same lambda prior for every quantile (broadcast scalar to length-K vector;
# tweak individual entries later if needed).
a_vec <- rep(a_lam, K)
b_vec <- rep(b_lam, K)

# ---------- Sampler hyper-parameters ----------
eps              <- 0.01
delta_vec_burn   <- rep(0.005, K) #0.001 - Normal #0.005 - Beta #0.01 - Mixed
# PC upper c(rep(0.010, 2), rep(0.009, K-10), rep(0.009, 2), 0.009, 0.009, 0.008, 0.007, 0.007, 0.007)
delta_vec_post   <- delta_vec_burn
gamma            <- 0.1
rho_admm         <- 1
max_iter_admm    <- 500
tol_admm         <- 1e-4

# Scale prior  sigma ~ IG(a_sigma, b_sigma)  (identical to individual sampler)
a_sigma <- 0.1
b_sigma <- 0.1

print_gibbs_every <- 500   # report every 500 iterations
print_admm_every  <- 0

gibbs_iter <- 6000
gibbs_burn <- 2000

# Quantiles to plot (subset of tau_vec); easy to edit
print_quantiles <- c(0.02)

# ---------- Synthetic data: noise switch identical to individual sampler ----------
if (noise_name == "Mixed") {
  y <- baseline + sapply(xs, function(xi) {
    if (runif(1) < xi) rnorm(1, -0.2, sqrt(0.5)) else rnorm(1, 0.2, sqrt(0.5))
  })
  Q_true <- matrix(NA, nrow = n_data, ncol = K)
  for (k in 1:K) {
    Q_true[, k] <- baseline + sapply(xs, function(xi) {
      cdf <- function(t) xi * pnorm(t, -0.2, sqrt(0.5)) +
                         (1 - xi) * pnorm(t,  0.2, sqrt(0.5))
      uniroot(function(t) cdf(t) - tau_vec[k], interval = c(-10, 10))$root
    })
  }
} else if (noise_name == "Beta") {
  y <- baseline + sapply(xs, function(xi) rbeta(1, 1, 11 - 10 * xi))
  Q_true <- matrix(NA, nrow = n_data, ncol = K)
  for (k in 1:K) {
    Q_true[, k] <- baseline + sapply(xs, function(xi) qbeta(tau_vec[k], 1, 11 - 10 * xi))
  }
} else {
  sd_noise <- (1 + xs^2) / 4
  y        <- baseline + rnorm(n_data, 0, sd_noise)
  Q_true <- matrix(NA, nrow = n_data, ncol = K)
  for (k in 1:K) {
    Q_true[, k] <- baseline + qnorm(tau_vec[k]) * sd_noise
  }
}

df <- tibble(x = x, y = y, baseline = baseline)
ggplot(df, aes(x)) +
  geom_point(aes(y = y), alpha = 0.45, size = 1, color = "blue") +
  geom_line(aes(y = baseline), linewidth = 1) +
  labs(title = "Data (points) and True Signal (line)",
       x = "x", y = "value") +
  theme_minimal()

# ---------- Helpers ----------
D_matrix <- function(n, order) {
  D <- diag(n)
  for (i in 1:(order + 1)) D <- diff(D)
  D
}

D <- D_matrix(n_data, order = order)
m <- nrow(D)

rho_tau <- function(u, tau) ifelse(u >= 0, tau * u, (tau - 1) * u)
true_check_loss <- function(u, tau) rho_tau(u, tau)

# Sigma-aware gradient/posterior, same forms as individual sampler
grad_smooth_check <- function(y, theta, tau, sigma, eps) {
  u <- (y - theta) / sigma
  sigmoid <- 1 / (1 + exp(u / eps))
  -(tau - sigmoid) / sigma
}

neg_log_posterior <- function(theta, y, tau, lambda, sigma, D, eps) {
  resid   <- y - theta
  f_theta <- sum((1/sigma) * true_check_loss(resid, tau))
  g_theta <- lambda * sum(abs(D %*% theta))
  f_theta + g_theta
}

# Inverse-Gamma sampler (X ~ Gamma(a, b)  =>  1/X ~ IG(a, b))
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

# ---------- ADMM wrapper (unchanged) ----------
admm_sequential <- function(w, D, gamma, lambda, lower_bound, upper_bound,
                            rho = 100, max_iter = 300, tol = 1e-4,
                            theta_init = NULL, v_init = NULL, w_init = NULL,
                            u1_init = NULL, u2_init = NULL,
                            verbose = FALSE, print_every = 0) {

  result <- prox_admm_optimized_fixed(
    w = w, lambda = lambda, gamma = gamma,
    lower_bound = lower_bound, upper_bound = upper_bound,
    k = order, rho = rho, max_iter = max_iter, tol = tol,
    theta_init = theta_init, v_init = v_init, w_init_param = w_init,
    u1_init = u1_init, u2_init = u2_init,
    verbose = verbose, print_every = print_every
  )

  list(theta = result$theta, v = result$v, w = result$w,
       u1 = result$u1, u2 = result$u2,
       r_primal_1 = result$r_primal_1, r_primal_2 = result$r_primal_2,
       r_dual_1 = result$r_dual_1, r_dual_2 = result$r_dual_2,
       theta_L1 = result$theta_L1, z_L1 = result$z_L1,
       w_L1 = result$w_L1, w_theta_L1 = result$w_theta_L1,
       iterations = result$iterations)
}

# ---------- P-MALA single sample (sigma-aware) ----------
pmala_single_sequential <- function(y, tau, lambda, sigma, gamma, D, delta, eps,
                                    theta_current, lower_bound, upper_bound,
                                    order = 2, rho_admm = 100,
                                    admm_state = NULL,
                                    verbose_admm = FALSE, print_admm_every = 0) {

  n <- length(y)

  grad      <- grad_smooth_check(y, theta_current, tau, sigma, eps)
  grad_norm <- sum(abs(grad))
  z         <- theta_current - delta * grad + sqrt(2 * delta) * rnorm(n)

  admm_result <- admm_sequential(z, D, gamma, lambda, lower_bound, upper_bound,
                                 rho = rho_admm, max_iter = max_iter_admm,
                                 tol = tol_admm,
                                 theta_init = z,
                                 v_init = NULL,
                                 w_init = admm_state$w,
                                 u1_init = NULL, u2_init = NULL,
                                 verbose = verbose_admm,
                                 print_every = print_admm_every)
  theta_proposal <- admm_result$theta

  log_post_curr <- -neg_log_posterior(theta_current,  y, tau, lambda, sigma, D, eps)
  log_post_prop <- -neg_log_posterior(theta_proposal, y, tau, lambda, sigma, D, eps)
  log_alpha     <- min(0, log_post_prop - log_post_curr)

  accepted  <- (log(runif(1)) < log_alpha)
  theta_new <- if (accepted) theta_proposal else theta_current

  list(theta = theta_new, accepted = accepted,
       admm_state = admm_result, grad_norm = grad_norm)
}

# ---------- Joint Gibbs: theta_k, lambda_k, sigma_k for k=1..K ----------
gibbs_sequential <- function(y, tau_vec, a_vec, b_vec,
                             a_sigma, b_sigma,
                             gamma, D, delta_vec_burn, delta_vec_post,
                             gibbs_iter, gibbs_burn, eps, order = 2,
                             rho_admm = 100, use_warm_start = TRUE,
                             print_every = 500,
                             verbose_admm = FALSE, print_admm_every = 0) {

  n <- length(y)
  K <- length(tau_vec)

  stopifnot(length(delta_vec_burn) == K)
  stopifnot(length(delta_vec_post) == K)
  stopifnot(length(a_vec) == K, length(b_vec) == K)

  theta_chain      <- array(NA, dim = c(gibbs_iter, n, K))
  lambda_chain     <- matrix(NA, nrow = gibbs_iter, ncol = K)
  sigma_chain      <- matrix(NA, nrow = gibbs_iter, ncol = K)
  acceptance_chain <- matrix(NA, nrow = gibbs_iter, ncol = K)
  grad_norms       <- matrix(NA, nrow = gibbs_iter, ncol = K)

  theta       <- matrix(0, nrow = n, ncol = K)
  lambda_vec  <- rep(10000, K)
  sigma_vec   <- rep(1,     K)
  admm_states <- vector("list", K)

  for (k in 1:K) theta[, k] <- y

  cat(sprintf("Iterations: %d (burn-in: %d)\n", gibbs_iter, gibbs_burn))
  cat(sprintf("Quantiles: %s\n", paste(sprintf("%.2f", tau_vec), collapse = ", ")))
  cat(sprintf("ADMM: rho=%.1f, max_iter=%d, tol=%.1e\n", rho_admm, max_iter_admm, tol_admm))
  cat(sprintf("Lambda prior (per quantile, all equal): a=%.4g, b=%.4g\n", a_vec[1], b_vec[1]))
  cat(sprintf("Sigma prior:  a_sigma=%.4g, b_sigma=%.4g\n", a_sigma, b_sigma))
  cat(sprintf("Reporting every %d iterations\n\n", print_every))

  start_time <- Sys.time()

  for (s in 1:gibbs_iter) {

    delta_vec_current <- if (s <= gibbs_burn) delta_vec_burn else delta_vec_post

    # ---- 1. theta_k | rest, sequentially with monotonicity box ----
    for (k in 1:K) {
      if (k == 1) {
        lower_bound <- rep(-Inf, n);    upper_bound <- theta[, k + 1]
      } else if (k == K) {
        lower_bound <- theta[, k - 1];  upper_bound <- rep(Inf, n)
      } else {
        lower_bound <- theta[, k - 1];  upper_bound <- theta[, k + 1]
      }

      pmala_result <- pmala_single_sequential(
        y, tau_vec[k], lambda_vec[k], sigma_vec[k],
        gamma, D, delta_vec_current[k], eps,
        theta[, k], lower_bound, upper_bound,
        order = order, rho_admm = rho_admm,
        admm_state = admm_states[[k]],
        verbose_admm = verbose_admm, print_admm_every = print_admm_every
      )

      theta[, k]             <- pmala_result$theta
      acceptance_chain[s, k] <- pmala_result$accepted
      grad_norms[s, k]       <- pmala_result$grad_norm
      admm_states[[k]]       <- pmala_result$admm_state
    }

    # ---- 2. lambda_k | theta_k ~ Gamma(a_k + m, b_k + ||D theta_k||_1) ----
    for (k in 1:K) {
      l1_norm        <- sum(abs(D %*% theta[, k]))
      lambda_vec[k]  <- rgamma(1, shape = a_vec[k] + m, rate = b_vec[k] + l1_norm)
    }

    # ---- 3. sigma_k | theta_k ~ IG(n + a_sigma, S_rho + b_sigma) ----
    for (k in 1:K) {
      S_rho        <- sum(true_check_loss(y - theta[, k], tau_vec[k]))
      sigma_vec[k] <- max(0.001, rinvgamma(1, shape = n + a_sigma, rate = S_rho + b_sigma))
    }

    theta_chain[s, , ] <- theta
    lambda_chain[s, ]  <- lambda_vec
    sigma_chain[s, ]   <- sigma_vec

    # ---- Reporting ----
    if (s %% print_every == 0 || s == 1) {
      elapsed     <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      iter_per_sec <- s / elapsed
      remaining   <- (gibbs_iter - s) / iter_per_sec

      window_size  <- min(print_every, s)
      window_start <- max(1, s - window_size + 1)
      accept_rates <- colMeans(acceptance_chain[window_start:s, , drop = FALSE], na.rm = TRUE)
      dtheta_l1    <- sapply(1:K, function(k) sum(abs(D %*% theta[, k])))

      cat(sprintf("\n[Gibbs Iter %d/%d]\n", s, gibbs_iter))
      cat("------------------------------------------------------------------\n")
      cat(sprintf("  Lambdas:       %s\n", paste(sprintf("%.2f", lambda_vec),  collapse = ", ")))
      cat(sprintf("  Sigmas:        %s\n", paste(sprintf("%.3f", sigma_vec),   collapse = ", ")))
      cat(sprintf("  ||D theta||_1: %s\n", paste(sprintf("%.2f", dtheta_l1),   collapse = ", ")))
      cat(sprintf("  Accept Rates:  %s\n", paste(sprintf("%.3f", accept_rates),collapse = ", ")))
      cat(sprintf("  Speed: %.2f it/s | ETA: %.0fs (%.1f min)\n",
                  iter_per_sec, remaining, remaining/60))
      cat("------------------------------------------------------------------\n")
    }
  }

  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("\nTotal time: %.1f seconds (%.2f iterations/sec)\n\n",
              total_time, gibbs_iter / total_time))

  # Post burn-in samples
  idx_post       <- (gibbs_burn + 1):gibbs_iter
  theta_samples  <- theta_chain[idx_post, , , drop = FALSE]
  lambda_samples <- lambda_chain[idx_post, , drop = FALSE]
  sigma_samples  <- sigma_chain[idx_post, , drop = FALSE]

  theta_mean   <- apply(theta_samples, c(2, 3), mean)
  theta_median <- apply(theta_samples, c(2, 3), median)

  # 95% credible bands (for MCIW)
  theta_lower <- apply(theta_samples, c(2, 3), quantile, probs = 0.025)
  theta_upper <- apply(theta_samples, c(2, 3), quantile, probs = 0.975)

  lambda_mean   <- colMeans(lambda_samples)
  lambda_median <- apply(lambda_samples, 2, median)
  sigma_mean    <- colMeans(sigma_samples)
  sigma_median  <- apply(sigma_samples, 2, median)

  acceptance_rates <- colMeans(acceptance_chain[idx_post, , drop = FALSE])

  list(
    theta_samples    = theta_samples,
    lambda_samples   = lambda_samples,
    sigma_samples    = sigma_samples,
    theta_chain      = theta_chain,
    lambda_chain     = lambda_chain,
    sigma_chain      = sigma_chain,
    theta_mean       = theta_mean,
    theta_median     = theta_median,
    theta_lower      = theta_lower,
    theta_upper      = theta_upper,
    lambda_mean      = lambda_mean,
    lambda_median    = lambda_median,
    sigma_mean       = sigma_mean,
    sigma_median     = sigma_median,
    acceptance_rates = acceptance_rates,
    acceptance_chain = acceptance_chain,
    grad_norms       = grad_norms
  )
}

# ---------- Run ----------
result <- gibbs_sequential(
  y = y, tau_vec = tau_vec,
  a_vec = a_vec, b_vec = b_vec,
  a_sigma = a_sigma, b_sigma = b_sigma,
  gamma = gamma, D = D,
  delta_vec_burn = delta_vec_burn,
  delta_vec_post = delta_vec_post,
  gibbs_iter = gibbs_iter, gibbs_burn = gibbs_burn,
  eps = eps, order = order,
  rho_admm = rho_admm, use_warm_start = TRUE,
  print_every = print_gibbs_every,
  verbose_admm = FALSE,
  print_admm_every = print_admm_every
)

# ---------- Diagnostics: lambda / sigma summaries ----------
cat("\nLambda estimates:\n")
for (k in 1:K) {
  cat(sprintf("  tau = %.3f: mean = %.3f, median = %.3f\n",
              tau_vec[k], result$lambda_mean[k], result$lambda_median[k]))
}

cat("\nSigma estimates:\n")
for (k in 1:K) {
  cat(sprintf("  tau = %.3f: mean = %.4f, median = %.4f\n",
              tau_vec[k], result$sigma_mean[k], result$sigma_median[k]))
}

cat("\nAcceptance rates:\n")
for (k in 1:K) {
  cat(sprintf("  tau = %.3f: %.3f\n", tau_vec[k], result$acceptance_rates[k]))
}
cat("\n")

# ---------- Output directories ----------
dir.create("plots",   showWarnings = FALSE)
dir.create("lambdas", showWarnings = FALSE)
dir.create("thetas",  showWarnings = FALSE)
dir.create("sigmas",  showWarnings = FALSE)

# ---------- Customizable subset plot: data + true + mean + median ----------
cat("Creating customizable subset plot for selected quantiles...\n")

# Match each requested tau to closest available tau_vec entry
match_idx <- sapply(print_quantiles, function(q) which.min(abs(tau_vec - q)))
match_idx <- unique(match_idx)
sub_taus  <- tau_vec[match_idx]

cat(sprintf("Plotting quantiles: %s\n",
            paste(sprintf("%.3f", sub_taus), collapse = ", ")))

df_sub <- do.call(rbind, lapply(seq_along(match_idx), function(j) {
  k <- match_idx[j]
  data.frame(
    x        = rep(x, 3),
    value    = c(Q_true[, k], result$theta_mean[, k], result$theta_median[, k]),
    type     = rep(c("True", "Mean", "Median"), each = n_data),
    quantile = factor(sprintf("tau = %.3f", tau_vec[k]),
                      levels = sprintf("tau = %.3f", sub_taus))
  )
}))

p_subset <- ggplot() +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y),
             color = "black", alpha = 0.45, size = 0.8) +
  geom_line(data = df_sub,
            aes(x = x, y = value, color = quantile, linetype = type),
            linewidth = 0.9) +
  scale_linetype_manual(values = c("True" = "solid",
                                   "Mean" = "dashed",
                                   "Median" = "dotted")) +
  labs(title = "Selected Quantile Estimates (Mean & Median) vs True",
       subtitle = paste("model =", model, " | noise =", noise_name),
       x = "x", y = "Quantile value",
       color = "Quantile", linetype = "Curve") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave("plots/subset_estimates.png", p_subset,
       width = 12, height = 7, dpi = 300)
print(p_subset)
cat("Exported: plots/subset_estimates.png\n\n")

# ---------- Final metrics table: MSE / MAD / MCIW for mean & median ----------
metrics_tbl <- data.frame(
  tau         = tau_vec,
  MSE_mean    = NA_real_,
  MAD_mean    = NA_real_,
  MCIW_mean   = NA_real_,
  MSE_median  = NA_real_,
  MAD_median  = NA_real_,
  MCIW_median = NA_real_
)

for (k in 1:K) {
  qt <- Q_true[, k]
  mean_est   <- result$theta_mean[, k]
  median_est <- result$theta_median[, k]
  width_k    <- mean(result$theta_upper[, k] - result$theta_lower[, k])

  metrics_tbl$MSE_mean[k]    <- mean((qt - mean_est)^2)
  metrics_tbl$MAD_mean[k]    <- mean(abs(qt - mean_est))
  metrics_tbl$MCIW_mean[k]   <- width_k     # CI is the same regardless of point estimator

  metrics_tbl$MSE_median[k]  <- mean((qt - median_est)^2)
  metrics_tbl$MAD_median[k]  <- mean(abs(qt - median_est))
  metrics_tbl$MCIW_median[k] <- width_k
}

cat("\n========================================================================\n")
cat(sprintf("FINAL METRICS  (model = %d, noise = %s)\n", model, noise_name))
cat("MCIW = mean width of 95%% credible interval (2.5%% - 97.5%%) across positions\n")
cat("========================================================================\n")
print(format(metrics_tbl, digits = 5, nsmall = 5), row.names = FALSE)
cat("========================================================================\n\n")
