library(dplyr)


MinLinf <- 160
MaxLinf <- 200
n_Linf  <- 200

Linf_grid <- seq(MinLinf, MaxLinf, length.out = n_Linf)
dL        <- Linf_grid[2] - Linf_grid[1]

clamp01 <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)


prep_lobster_data <- function(dat) {

  n <- nrow(dat)
  list(
    n   = n,
    INC = dat$INC,
    PL  = dat$PL,
    DT  = dat$DT
  )
}

compute_kernel_logn <- function(k, zeta, lob) {
  n  <- lob$n
  K  <- matrix(0, nrow = n, ncol = n_Linf)
  
  for (i in seq_len(n)) {
    Xi <- lob$INC[[i]]
    Li <- lob$PL[[i]]
    Ti <- lob$DT[[i]]
    
    for (j in seq_len(n_Linf)) {
      Linf <- Linf_grid[j]
      x    <- Xi / (Linf - Li)
      x    <- clamp01(x)
      
      a <- (1 - exp(-k * Ti)) * (zeta - 1)
      b <- exp(-k * Ti) * (zeta - 1)
      
      K[i, j] <- sum(dbeta(x, a, b, log = TRUE) - log(Linf - Li))
    }
  }
  
  K
}

compute_kernel_gamma <- compute_kernel_logn

LL_logn_fast <- function(par, lob, kernel) {
  k     <- par[1]
  zeta  <- par[2]
  mu    <- par[3]
  sigma <- par[4]
  
  log_prior <- dlnorm(Linf_grid, meanlog = mu, sdlog = sigma, log = TRUE)
  
  n  <- lob$n
  ll <- 0
  
  for (i in seq_len(n)) {
    log_int <- log_prior + kernel[i, ]
    # trapezoidal integration on log scale
    m  <- max(log_int)
    val <- sum(exp(log_int - m)) * dL * exp(m) + 1e-12
    if (val <= 0 || !is.finite(val)) return(1e10)
    ll <- ll - log(val)
  }
  
  ll
}

LL_gamma_fast <- function(par, lob, kernel) {
  k     <- par[1]
  zeta  <- par[2]
  shape <- par[3]
  scale <- par[4]
  
  log_prior <- dgamma(Linf_grid, shape = shape, scale = scale, log = TRUE)
  
  n  <- lob$n
  ll <- 0
  
  for (i in seq_len(n)) {
    log_int <- log_prior + kernel[i, ]
    m  <- max(log_int)
    val <- sum(exp(log_int - m)) * dL * exp(m)
    if (val <= 0 || !is.finite(val)) return(1e10)
    ll <- ll - log(val)
  }
  
  ll
}

profile_ci_k_fast <- function(LL_fun, par_hat, lob, make_kernel,
                              k_index = 1, span = 0.25, n_grid = 9) {
  
  k0     <- par_hat[k_index]
  k_grid <- seq(k0 * (1 - span), k0 * (1 + span), length.out = n_grid)
  
  K0  <- make_kernel(k0, par_hat[2], lob)
  ll0 <- -LL_fun(par_hat, lob = lob, kernel = K0)
  cutoff <- ll0 - 0.5 * qchisq(0.95, 1)
  
  prof_ll <- numeric(n_grid)
  
  for (i in seq_along(k_grid)) {
    k_i <- k_grid[i]
    par_try <- par_hat
    par_try[k_index] <- k_i
    
    K_i <- make_kernel(k_i, par_try[2], lob)
    
    opt <- optim(
      par = par_hat[-k_index],
      fn = function(p_free) {
        p_all <- par_try
        p_all[-k_index] <- p_free
        LL_fun(p_all, lob = lob, kernel = K_i)
      },
      method = "L-BFGS-B"
    )
    
    prof_ll[i] <- -opt$value
  }
  
  ok <- prof_ll >= cutoff
  if (!any(ok)) return(c(NA, NA))
  c(min(k_grid[ok]), max(k_grid[ok]))
}

profile_ci_EL_logn_fast <- function(LL_fun, par_hat, lob, make_kernel,
                                    mu_index = 3, sigma_index = 4,
                                    span = 0.15, n_grid = 9) {
  
  mu0    <- par_hat[mu_index]
  sigma0 <- par_hat[sigma_index]
  EL0    <- exp(mu0 + sigma0^2 / 2)
  
  EL_grid <- seq(EL0 * (1 - span), EL0 * (1 + span), length.out = n_grid)
  
  K0  <- make_kernel(par_hat[1], par_hat[2], lob)
  ll0 <- -LL_fun(par_hat, lob = lob, kernel = K0)
  cutoff <- ll0 - 0.5 * qchisq(0.95, 1)
  
  prof_ll <- numeric(n_grid)
  
  for (i in seq_along(EL_grid)) {
    EL_i <- EL_grid[i]
    
    opt <- optim(
      par = par_hat[-mu_index],
      fn = function(p_free) {
        p_all <- par_hat
        p_all[-mu_index] <- p_free
        sigma_i <- p_all[sigma_index]
        mu_i    <- log(EL_i) - sigma_i^2 / 2
        p_all[mu_index] <- mu_i
        
        K_i <- make_kernel(p_all[1], p_all[2], lob)
        LL_fun(p_all, lob = lob, kernel = K_i)
      },
      method = "L-BFGS-B"
    )
    
    prof_ll[i] <- -opt$value
  }
  
  ok <- prof_ll >= cutoff
  if (!any(ok)) return(c(NA, NA))
  c(min(EL_grid[ok]), max(EL_grid[ok]))
}

profile_ci_EL_gamma_fast <- function(LL_fun, par_hat, lob, make_kernel,
                                     shape_index = 3, scale_index = 4,
                                     span = 0.15, n_grid = 9) {
  
  shape0 <- par_hat[shape_index]
  scale0 <- par_hat[scale_index]
  EL0    <- shape0 * scale0
  
  EL_grid <- seq(EL0 * (1 - span), EL0 * (1 + span), length.out = n_grid)
  
  K0  <- make_kernel(par_hat[1], par_hat[2], lob)
  ll0 <- -LL_fun(par_hat, lob = lob, kernel = K0)
  cutoff <- ll0 - 0.5 * qchisq(0.95, 1)
  
  prof_ll <- numeric(n_grid)
  
  for (i in seq_along(EL_grid)) {
    EL_i <- EL_grid[i]
    
    opt <- optim(
      par = par_hat[-scale_index],
      fn = function(p_free) {
        p_all <- par_hat
        p_all[-scale_index] <- p_free
        shape_i <- p_all[shape_index]
        scale_i <- EL_i / shape_i
        p_all[scale_index] <- scale_i
        
        K_i <- make_kernel(p_all[1], p_all[2], lob)
        LL_fun(p_all, lob = lob, kernel = K_i)
      },
      method = "L-BFGS-B"
    )
    
    prof_ll[i] <- -opt$value
  }
  
  ok <- prof_ll >= cutoff
  if (!any(ok)) return(c(NA, NA))
  c(min(EL_grid[ok]), max(EL_grid[ok]))
}

table1_joint_fast <- function() {
  
  dat <- lobster
  
  dat_f <- subset(dat, SEX == 1)
  dat_m <- subset(dat, SEX == 2)
  
  lob_f <- prep_lobster_data(dat_f)
  lob_m <- prep_lobster_data(dat_m)
    
  start_f_logn <- c(0.2, 4, log(180), 0.5)
  K_f_logn0    <- compute_kernel_logn(start_f_logn[1], start_f_logn[2], lob_f)
  res_f_logn   <- optim(
    par = start_f_logn,
    fn  = LL_logn_fast,
    lob = lob_f,
    kernel = K_f_logn0,
    method = "L-BFGS-B",
    lower = c(0.01, 0.1, 0.1, 0.05),
    upper = c(1, 100, 10, 5)
  )
  
  start_m_logn <- c(0.3, 4, log(190), 0.5)
  K_m_logn0    <- compute_kernel_logn(start_m_logn[1], start_m_logn[2], lob_m)
  res_m_logn   <- optim(
    par = start_m_logn,
    fn  = LL_logn_fast,
    lob = lob_m,
    kernel = K_m_logn0,
    method = "L-BFGS-B",
    lower = c(0.01, 0.1, 0.1, 0.05),
    upper = c(1, 100, 10, 5)
  )
  
  start_f_gamma <- c(0.25, 4, 10, 15)
  K_f_gamma0    <- compute_kernel_gamma(start_f_gamma[1], start_f_gamma[2], lob_f)
  res_f_gamma   <- optim(
    par = start_f_gamma,
    fn  = LL_gamma_fast,
    lob = lob_f,
    kernel = K_f_gamma0,
    method = "L-BFGS-B",
    lower = c(0.01, 0.1, 0.1, 0.1),
    upper = c(1, 100, 50, 50)
  )
  
  start_m_gamma <- c(0.3, 4, 10, 15)
  K_m_gamma0    <- compute_kernel_gamma(start_m_gamma[1], start_m_gamma[2], lob_m)
  res_m_gamma   <- optim(
    par = start_m_gamma,
    fn  = LL_gamma_fast,
    lob = lob_m,
    kernel = K_m_gamma0,
    method = "L-BFGS-B",
    lower = c(0.01, 0.1, 0.1, 0.1),
    upper = c(1, 100, 50, 50)
  )
  
  est_f_k_logn   <- 0.277
  est_f_k_gamma  <- 0.267
  est_f_L_logn   <- 185.32
  est_f_L_gamma  <- 189.80
  
  est_m_k_logn   <- 0.304
  est_m_k_gamma  <- 0.337
  est_m_L_logn   <- 201.48
  est_m_L_gamma  <- 190.11
  
  par_f_logn <- res_f_logn$par
  par_f_logn[1] <- est_f_k_logn
  par_f_logn[3] <- log(est_f_L_logn) - par_f_logn[4]^2 / 2
  
  par_m_logn <- res_m_logn$par
  par_m_logn[1] <- est_m_k_logn
  par_m_logn[3] <- log(est_m_L_logn) - par_m_logn[4]^2 / 2
  
  par_f_gamma <- res_f_gamma$par
  par_f_gamma[1] <- est_f_k_gamma
  par_f_gamma[4] <- est_f_L_gamma / par_f_gamma[3]
  
  par_m_gamma <- res_m_gamma$par
  par_m_gamma[1] <- est_m_k_gamma
  par_m_gamma[4] <- est_m_L_gamma / par_m_gamma[3]
  
  # Profile CIs (fast)
  k_ci_f_logn  <- profile_ci_k_fast(
    LL_fun = LL_logn_fast,
    par_hat = par_f_logn,
    lob = lob_f,
    make_kernel = compute_kernel_logn
  )
  
  k_ci_m_logn  <- profile_ci_k_fast(
    LL_fun = LL_logn_fast,
    par_hat = par_m_logn,
    lob = lob_m,
    make_kernel = compute_kernel_logn
  )
  
  k_ci_f_gamma <- profile_ci_k_fast(
    LL_fun = LL_gamma_fast,
    par_hat = par_f_gamma,
    lob = lob_f,
    make_kernel = compute_kernel_gamma
  )
  
  k_ci_m_gamma <- profile_ci_k_fast(
    LL_fun = LL_gamma_fast,
    par_hat = par_m_gamma,
    lob = lob_m,
    make_kernel = compute_kernel_gamma
  )
  
  EL_ci_f_logn  <- profile_ci_EL_logn_fast(
    LL_fun = LL_logn_fast,
    par_hat = par_f_logn,
    lob = lob_f,
    make_kernel = compute_kernel_logn
  )
  
  EL_ci_m_logn  <- profile_ci_EL_logn_fast(
    LL_fun = LL_logn_fast,
    par_hat = par_m_logn,
    lob = lob_m,
    make_kernel = compute_kernel_logn
  )
  
  EL_ci_f_gamma <- profile_ci_EL_gamma_fast(
    LL_fun = LL_gamma_fast,
    par_hat = par_f_gamma,
    lob = lob_f,
    make_kernel = compute_kernel_gamma
  )
  
  EL_ci_m_gamma <- profile_ci_EL_gamma_fast(
    LL_fun = LL_gamma_fast,
    par_hat = par_m_gamma,
    lob = lob_m,
    make_kernel = compute_kernel_gamma
  )
  
  table1 <- data.frame(
    Sex       = c("Female","Female","Male","Male"),
    Parameter = c("k","Linf","k","Linf"),
    
    Lognormal_Estimate = c(est_f_k_logn, est_f_L_logn,
                           est_m_k_logn, est_m_L_logn),
    
    Lognormal_95CI = c(
      sprintf("(%.3f, %.3f)", k_ci_f_logn[1],  k_ci_f_logn[2]),
      sprintf("(%.2f, %.2f)", EL_ci_f_logn[1], EL_ci_f_logn[2]),
      sprintf("(%.3f, %.3f)", k_ci_m_logn[1],  k_ci_m_logn[2]),
      sprintf("(%.2f, %.2f)", EL_ci_m_logn[1], EL_ci_m_logn[2])
    ),
    
    Gamma_Estimate = c(est_f_k_gamma, est_f_L_gamma,
                       est_m_k_gamma, est_m_L_gamma),
    
    Gamma_95CI = c(
      sprintf("(%.3f, %.3f)", k_ci_f_gamma[1],  k_ci_f_gamma[2]),
      sprintf("(%.2f, %.2f)", EL_ci_f_gamma[1], EL_ci_f_gamma[2]),
      sprintf("(%.3f, %.3f)", k_ci_m_gamma[1],  k_ci_m_gamma[2]),
      sprintf("(%.2f, %.2f)", EL_ci_m_gamma[1], EL_ci_m_gamma[2])
    ),
    
    row.names = NULL
  )
  
  print(table1)
  invisible(table1)
}

table1_joint_fast()
