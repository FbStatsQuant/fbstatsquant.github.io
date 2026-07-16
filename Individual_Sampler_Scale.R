library(glmgen)
library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)

set.seed(800)
n_data <- 400
x <- seq(1, n_data, 1)


x_idx  <- x
xs <- x / n_data


# Parameters
tau <- 0.5
order <- 2
eps <- 0.01
delta <- 0.04

#Beta   0.005
#Normal 0.01
#Mixed 0.04
gamma <- 0.1

gibbs_iter <-6000
gibbs_burn <- 1000


model <- 3
noise_name <- "Mixed"


if (model == 1){
  m1 <- n_data/100
  baseline <- 2.5 * (x_idx >= 1   & x_idx <= m1*20)  +
    1.0 * (x_idx >= m1*20+1  & x_idx <= m1*40) +
    3.5 * (x_idx >= m1*40+1 & x_idx <= m1*60) +
    1.5 * (x_idx >= m1*60+1 & x_idx <= m1*100)
  C_1 <-0.5
  a <- C_1*n_data
  b <- n_data
  order <- 0
  
  a <- 1
  b <- C_1*n_data^(2/3)

} else if (model == 2){
  c_1 <- 1
  baseline <- 2 + sin(4 * xs - 2) + 2 * exp(-30 * (4 * xs - 2)^2)
  a <- 1
  b <- (1/c_1)*sqrt(n_data)
  order <- 1
  
  #a <- C_1*n_data
  #b <- n_data^(3/5)
  
  #a <- 1
  #b <- 1
  #order <- 2
} else if (model == 3){
  c_1 <- 1
  baseline <- sin(2 * pi * x / n_data)
  a <- 1
  b <- 1
  order <- 2
  
  #a <- 1
  #b <- 1
  #order <- 2
}

sd_noise <- (1 + xs^2) / 4

if (noise_name == "Mixed") {
  y <- baseline + sapply(xs, function(xi) {
    if (runif(1) < xi) rnorm(1, -0.2, sqrt(0.5)) else rnorm(1, 0.2, sqrt(0.5))
  })
  Q_true <- baseline + sapply(xs, function(xi) {
    cdf <- function(t) xi * pnorm(t, -0.2, sqrt(0.5)) + (1 - xi) * pnorm(t, 0.2, sqrt(0.5))
    uniroot(function(t) cdf(t) - tau, interval = c(-10, 10))$root
  })
} else if (noise_name == "Beta") {
  y      <- baseline + sapply(xs, function(xi) rbeta(1, 1, 11 - 10 * xi))
  Q_true <- baseline + sapply(xs, function(xi) qbeta(tau, 1, 11 - 10 * xi))
} else {
  sd_noise <- (1 + xs^2) / 4
  y        <- baseline + rnorm(n_data, 0, sd_noise)
  Q_true   <- baseline + qnorm(tau) * sd_noise
}

df <- tibble(x = x, y = y, baseline = baseline)
initial_plot <- ggplot(df, aes(x)) +
  geom_point(aes(y = y), alpha = 0.45, size = 1, color = "blue") +
  geom_line(aes(y = baseline), linewidth = 1) +
  labs(title = "Data (points) and True Signal (line)", x = "x", y = "value") +
  theme_minimal()

#windows()
#print(initial_plot)

# IG hyperparameters for sigma
a_sigma <- 0.01
b_sigma <- 0.01

cat("Mean lambda prior: ", a/b, "\n")
cat("Sd lambda prior: ", sqrt(a/b^2), "\n")

cat("=== THREE-BLOCK GIBBS SAMPLER: theta, lambda, sigma ===\n")
cat("Inner P-MALA: 1 sample per Gibbs iteration (no burn-in)\n")
cat("Outer Gibbs: tracking acceptance at Gibbs level\n")
cat("Gamma prior: lambda ~ Gamma(", a, ",", b, ")\n")
cat("IG prior: sigma ~ IG(", a_sigma, ",", b_sigma, ")\n\n")

# FUNCTIONS
D_matrix <- function(n, order) {
  D <- diag(n)
  for (i in 1:(order + 1)) {
    D <- diff(D)
  }
  D
}

D <- D_matrix(n_data, order = order)

rho_tau <- function(u, tau) ifelse(u >= 0, tau*u, (tau-1)*u)

true_check_loss <- function(u, tau) {
  loss <- rho_tau(u, tau)
  (loss)
}

smooth_check_loss <- function(u, tau, eps = eps) {
  tau * u + eps * log(1 + exp(-u / eps))
}

# UPDATED: sigma is now an explicit argument
grad_smooth_check <- function(y, theta, tau, sigma, eps = eps) {
  u <- (y - theta) / sigma
  sigmoid <- 1 / (1 + exp(u / eps))
  -(tau - sigmoid) / sigma
}

# UPDATED: sigma is now an explicit argument
neg_log_posterior_true <- function(theta, y, tau, lambda, sigma, D, eps) {
  resid <- y - theta
  f_theta <- sum((1/sigma) * true_check_loss(resid, tau))
  g_theta <- lambda * sum(abs(D %*% theta))
  f_theta + g_theta
}

prox_g <- function(z, D, gamma, lambda, order = 2) {
  eff <- gamma * lambda
  fit <- glmgen::trendfilter(x = seq_along(z), y = z, k = order, lambda = eff,
                             control = trendfilter.control.list(max_iter = 800L))
  as.numeric(fit$beta)
}

grad_f_plus_g <- function(y, theta, tau, sigma, eps, gamma, D, lambda, order = 2) {
  grad_f_theta <- grad_smooth_check(y, theta, tau, sigma, eps)
  x_tilde <- theta - gamma * grad_f_theta
  prox_step <- prox_g(x_tilde, D, gamma, lambda, order = order)
  (theta - prox_step) / gamma
}

# P-MALA SINGLE SAMPLE with sigma
pmala_single_sample <- function(y, tau = 0.5, lambda = 50, sigma = 1, gamma = 0.1,
                                D, delta = 0.1, eps = 0.01,
                                theta_current = NULL, order = 2) {
  n <- length(y)
  
  if (is.null(theta_current)) {
    theta_current <- rep(0, n)
  }
  
  # Generate proposal
  grad <- grad_smooth_check(y, theta_current, tau, sigma, eps)
  z <- theta_current - delta * grad + sqrt(2 * delta) * rnorm(n)
  theta_proposal <- prox_g(z, D, gamma, lambda, order = order)
  
  # Compute acceptance probability
  log_posterior_current  <- -neg_log_posterior_true(theta_current,  y, tau, lambda, sigma, D, eps)
  log_posterior_proposal <- -neg_log_posterior_true(theta_proposal, y, tau, lambda, sigma, D, eps)
  log_posterior_ratio <- log_posterior_proposal - log_posterior_current
  
  grad_psi_current <- grad_f_plus_g(y, theta_current, tau, sigma, eps, gamma, D, lambda, order)
  mu_current <- theta_current - delta * grad_psi_current
  
  grad_psi_proposal <- grad_f_plus_g(y, theta_proposal, tau, sigma, eps, gamma, D, lambda, order)
  mu_proposal <- theta_proposal - delta * grad_psi_proposal
  
  log_q_proposal_given_current <- -sum((theta_proposal - mu_current)^2) / (4 * delta)
  log_q_current_given_proposal <- -sum((theta_current - mu_proposal)^2) / (4 * delta)
  log_proposal_ratio <- log_q_current_given_proposal - log_q_proposal_given_current
  
  log_alpha <- min(0, log_posterior_ratio + log_proposal_ratio)
  
  # Accept/reject
  accepted <- (log(runif(1)) < log_alpha)
  
  if (accepted) {
    theta_new <- theta_proposal
  } else {
    theta_new <- theta_current
  }
  
  list(theta = theta_new, accepted = accepted)
}

# Inverse-Gamma sampler: if X ~ Gamma(a, rate=b), then 1/X ~ IG(a, b)
rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

gibbs_three_block <- function(y, tau = 0.5,
                              a = 2, b = 0.002,
                              a_sigma = 0.1, b_sigma = 0.1,
                              gamma = 0.1, D, delta = 0.1,
                              gibbs_iter = 30, gibbs_burn = 10,
                              eps = 0.01, theta_init = NULL,
                              sigma_init = 1, lambda_init = 20,
                              order = 2) {
  
  n <- length(y)
  
  gibbs_post <- gibbs_iter - gibbs_burn
  theta_gibbs_samples  <- matrix(NA, nrow = gibbs_post, ncol = n)
  lambda_gibbs_samples <- numeric(gibbs_post)
  sigma_gibbs_samples  <- numeric(gibbs_post)
  
  theta_gibbs_chain  <- matrix(NA, nrow = gibbs_iter, ncol = n)
  lambda_gibbs_chain <- numeric(gibbs_iter)
  sigma_gibbs_chain  <- numeric(gibbs_iter)
  
  acceptance_count <- 0
  
  if (is.null(theta_init)) {
    theta_current <- y
  } else {
    theta_current <- theta_init
  }
  lambda_current <- lambda_init
  sigma_current  <- sigma_init
  
  cat("Starting three-block Gibbs sampler (theta, lambda, sigma)...\n")
  cat("Initial: lambda =", lambda_current, ", sigma =", sigma_current, "\n\n")
  
  m <- nrow(D)
  
  for (i in 1:gibbs_iter) {
    
    # --- 1. Sample theta | y, lambda, sigma via P-MALA ---
    pmala_result <- pmala_single_sample(
      y = y,
      tau = tau,
      lambda = lambda_current,
      sigma = sigma_current,
      gamma = gamma,
      D = D,
      delta = delta,
      eps = eps,
      theta_current = theta_current,
      order = order
    )
    
    theta_current <- pmala_result$theta
    
    if (pmala_result$accepted) {
      acceptance_count <- acceptance_count + 1
    }
    
    # --- 2. Sample lambda | y, theta, sigma ~ Gamma(a + m, b + S_D) ---
    S_D <- sum(abs(D %*% theta_current))
    lambda_current <- rgamma(1, shape = a + m, rate = b + S_D)
    
    # --- 3. Sample sigma | y, theta, lambda ~ IG(n + a_sigma, S_rho + b_sigma) ---
    S_rho <- sum(true_check_loss(y - theta_current, tau))
    sigma_current <- rinvgamma(1, shape = n + a_sigma, rate = S_rho + b_sigma)
    
    if (i %% 500 == 0) {
      acc_rate <- acceptance_count / i
      cat("Iter", i,
          "| lambda =", round(lambda_current, 3),
          "| sigma =", round(sigma_current, 4),
          "| S_D =", round(S_D, 4),
          "| S_rho =", round(S_rho, 2),
          "| Acc Rate:", round(acc_rate, 3), "\n")
    }
    
    # Store chains
    theta_gibbs_chain[i, ]  <- theta_current
    lambda_gibbs_chain[i]   <- lambda_current
    sigma_gibbs_chain[i]    <- sigma_current
    
    if (i > gibbs_burn) {
      sample_idx <- i - gibbs_burn
      theta_gibbs_samples[sample_idx, ]  <- theta_current
      lambda_gibbs_samples[sample_idx]   <- lambda_current
      sigma_gibbs_samples[sample_idx]    <- sigma_current
    }
  }
  
  compute_mode <- function(x) {
    d <- density(x)
    d$x[which.max(d$y)]
  }
  
  theta_mean_estimate   <- colMeans(theta_gibbs_samples)
  theta_median_estimate <- apply(theta_gibbs_samples, 2, median)
  theta_mode_estimate   <- apply(theta_gibbs_samples, 2, compute_mode)
  lambda_final_estimate <- mean(lambda_gibbs_samples)
  sigma_final_estimate  <- mean(sigma_gibbs_samples)
  
  list(
    theta_estimate        = theta_mean_estimate,
    theta_mean_estimate   = theta_mean_estimate,
    theta_median_estimate = theta_median_estimate,
    theta_mode_estimate   = theta_mode_estimate,
    lambda_estimate       = lambda_final_estimate,
    sigma_estimate        = sigma_final_estimate,
    theta_samples         = theta_gibbs_samples,
    lambda_samples        = lambda_gibbs_samples,
    sigma_samples         = sigma_gibbs_samples,
    theta_chain           = theta_gibbs_chain,
    lambda_chain          = lambda_gibbs_chain,
    sigma_chain           = sigma_gibbs_chain,
    acceptance_rate       = acceptance_count / gibbs_iter
  )
}


theta_init_smart <- y

result_gibbs <- gibbs_three_block(
  y = y,
  tau = tau,
  a = a, b = b,
  a_sigma = a_sigma, b_sigma = b_sigma,
  gamma = gamma,
  D = D,
  delta = delta,
  gibbs_iter = gibbs_iter,
  gibbs_burn = gibbs_burn,
  eps = eps,
  theta_init = theta_init_smart,
  sigma_init = 1,
  lambda_init = 20,
  order = order
)

cat("\n=== FINAL GIBBS RESULTS (three-block: theta, lambda, sigma) ===\n")
cat("Final lambda estimate (posterior mean):", round(result_gibbs$lambda_estimate, 4), "\n")
cat("Lambda posterior std:", round(sd(result_gibbs$lambda_samples), 4), "\n")
cat("Final sigma estimate (posterior mean):", round(result_gibbs$sigma_estimate, 4), "\n")
cat("Sigma posterior std:", round(sd(result_gibbs$sigma_samples), 4), "\n")
cat("Overall acceptance rate:", round(result_gibbs$acceptance_rate, 3), "\n")


mse_mean   <- mean((Q_true - result_gibbs$theta_mean_estimate)^2)
mse_median <- mean((Q_true - result_gibbs$theta_median_estimate)^2)
mse_mode   <- mean((Q_true - result_gibbs$theta_mode_estimate)^2)

mad_mean   <- mean(abs(Q_true - result_gibbs$theta_mean_estimate))
mad_median <- mean(abs(Q_true - result_gibbs$theta_median_estimate))
mad_mode   <- mean(abs(Q_true - result_gibbs$theta_mode_estimate))


df_plot <- tibble(
  x = x,
  y = y,
  theoretical = Q_true,
  mean_estimate   = result_gibbs$theta_mean_estimate,
  median_estimate = result_gibbs$theta_median_estimate,
  mode_estimate   = result_gibbs$theta_mode_estimate
)

# PLOT 1: MEAN ESTIMATE
p_mean <- ggplot(df_plot, aes(x = x)) +
  geom_point(aes(y = y), alpha = 0.3, color = "gray60", size = 0.5) +
  geom_line(aes(y = theoretical), color = "black", size = 1.2) +
  geom_line(aes(y = mean_estimate), color = "red", size = 1.0, linetype = "dashed") +
  labs(
    title = "Three-Block Gibbs (theta, lambda, sigma): MEAN Estimate",
    subtitle = paste("MSE:", round(mse_mean, 4),
                     "| MAD:", round(mad_mean, 4),
                     "| lambda:", round(result_gibbs$lambda_estimate, 2),
                     "| sigma:", round(result_gibbs$sigma_estimate, 3)),
    x = "x", y = "Quantile Value"
  ) +
  theme_minimal()

# PLOT 2: MEDIAN ESTIMATE
p_median <- ggplot(df_plot, aes(x = x)) +
  geom_point(aes(y = y), alpha = 0.3, color = "gray60", size = 0.5) +
  geom_line(aes(y = theoretical), color = "black", size = 1.2) +
  geom_line(aes(y = median_estimate), color = "blue", size = 1.0, linetype = "dashed") +
  labs(
    title = "Three-Block Gibbs (theta, lambda, sigma): MEDIAN Estimate",
    subtitle = paste("MSE:", round(mse_median, 4),
                     "| MAD:", round(mad_median, 4),
                     "| lambda:", round(result_gibbs$lambda_estimate, 2),
                     "| sigma:", round(result_gibbs$sigma_estimate, 3)),
    x = "x", y = "Quantile Value"
  ) +
  theme_minimal()

# PLOT 3: MODE ESTIMATE
p_mode <- ggplot(df_plot, aes(x = x)) +
  geom_point(aes(y = y), alpha = 0.3, color = "gray60", size = 0.5) +
  geom_line(aes(y = theoretical), color = "black", size = 1.2) +
  geom_line(aes(y = mode_estimate), color = "green", size = 1.0, linetype = "dashed") +
  labs(
    title = "Three-Block Gibbs (theta, lambda, sigma): MODE Estimate",
    subtitle = paste("MSE:", round(mse_mode, 4),
                     "| MAD:", round(mad_mode, 4),
                     "| lambda:", round(result_gibbs$lambda_estimate, 2),
                     "| sigma:", round(result_gibbs$sigma_estimate, 3)),
    x = "x", y = "Quantile Value"
  ) +
  theme_minimal()

# Trace plots for sigma and lambda (useful diagnostics)
df_trace <- tibble(
  iter   = 1:gibbs_iter,
  lambda = result_gibbs$lambda_chain,
  sigma  = result_gibbs$sigma_chain
)

p_lambda_trace <- ggplot(df_trace, aes(x = iter, y = lambda)) +
  geom_line(color = "darkred", alpha = 0.7) +
  geom_vline(xintercept = gibbs_burn, linetype = "dashed") +
  labs(title = "Lambda trace", x = "Iteration", y = "lambda") +
  theme_minimal()

p_sigma_trace <- ggplot(df_trace, aes(x = iter, y = sigma)) +
  geom_line(color = "darkblue", alpha = 0.7) +
  geom_vline(xintercept = gibbs_burn, linetype = "dashed") +
  labs(title = "Sigma trace", x = "Iteration", y = "sigma") +
  theme_minimal()

windows(); print(p_mean)
windows(); print(p_median)
windows(); print(p_mode)
#windows(); print(p_lambda_trace)
#windows(); print(p_sigma_trace)


cat("\nMSE and MAD (vs Q_true):\n")
cat("Mean estimate   MSE:", round(mse_mean,   6), "  MAD:", round(mad_mean,   6), "\n")
cat("Median estimate MSE:", round(mse_median, 6), "  MAD:", round(mad_median, 6), "\n")
cat("Mode estimate   MSE:", round(mse_mode,   6), "  MAD:", round(mad_mode,   6), "\n\n")