test_table2_bl <- function() {
  
  cat("Running test for table2 (BL model)...\n")
  
  dat <- lobster
  
  # 1. Log-posterior for BL model
  ff_ind <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
    
    if (!is.finite(Linf) || Linf <= max(Li) + 1e-8) return(-1e10)
    
    # Lognormal prior on Linf
    logPrior <- -log(Linf) - log(sigma) - 0.5 * log(2 * pi) -
      ((log(Linf) - mu)^2 / (2 * sigma^2))
    
    out <- logPrior
    
    for (j in seq_along(Xi)) {
      
      denom <- Linf - Li[j]
      if (denom <= 1e-8) return(-1e10)
      
      x <- Xi[j] / denom
      if (x <= 1e-8 || x >= 1 - 1e-8) return(-1e10)
      
      a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
      b <- exp(-k * Ti[j]) * (zeta - 1)
      if (a <= 1e-8 || b <= 1e-8) return(-1e10)
      
      out <- out + dbeta(x, a, b, log = TRUE) - log(denom)
    }
    
    out
  }
  
  kernel_norm <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
    exp(ff_ind(Linf, k, zeta, mu, sigma, Xi, Ti, Li) - Mi)
  }
  
  max_over_Linf <- function(fun, lower, upper) {
    vals <- c(
      optimize(fun, interval = c(lower, upper), maximum = TRUE)$objective,
      fun(lower), fun(upper), fun((lower + upper)/2)
    )
    max(vals)
  }
  
  safe_integrate <- function(fun, lower, upper, n = 300) {
    xs <- seq(lower, upper, length.out = n)
    vals <- sapply(xs, fun)
    dx <- (upper - lower) / (n - 1)
    sum(vals) * dx
  }
  
  # 2. Marginal log-likelihood
  LL <- function(theta, dat, MinLinf, MaxLinf) {
    
    k     <- theta[1]
    zeta  <- theta[2]
    mu    <- theta[3]
    sigma <- theta[4]
    
    if (k <= 0 || zeta <= 1.01 || sigma <= 0) return(1e12)
    
    LLtot <- 0
    
    for (id in unique(dat$LOBSTER)) {
      
      dati <- dat[dat$LOBSTER == id, ]
      Xi <- dati$INC
      Li <- dati$PL
      Ti <- dati$INT / 365.25
      
      Mi <- max_over_Linf(
        function(Ls) ff_ind(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      
      val <- safe_integrate(
        function(Ls) kernel_norm(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      
      if (!is.finite(val) || val <= 0) return(1e12)
      
      LLtot <- LLtot + log(val) + Mi
    }
    
    -LLtot
  }
  
  # 3. SEs for k and E[Linf] only (delta method)
  get_se_k_Linf <- function(opt, dat, MinLinf, MaxLinf) {
    
    hess <- optimHess(
      opt$par, fn = LL,
      dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf
    )
    
    vcov <- tryCatch(
      solve(hess),
      error = function(e) matrix(NA, 4, 4)
    )
    
    # SE for k
    se_k <- sqrt(vcov[1, 1])
    
    # Delta method for E[Linf] = exp(mu + sigma^2 / 2)
    mu    <- opt$par[3]
    sigma <- opt$par[4]
    
    ELinf <- exp(mu + 0.5 * sigma^2)
    
    grad <- c(
      d_mu    = ELinf,
      d_sigma = sigma * ELinf
    )
    
    var_ELinf <- grad %*% vcov[3:4, 3:4] %*% grad
    se_ELinf  <- sqrt(var_ELinf)
    
    list(se_k = se_k, se_Linf = se_ELinf)
  }
  
  # 4. FEMALE
  dat_f <- dat[dat$SEX == 1, ]
  MinLinf <- max(dat_f$PL + dat_f$INC) + 0.1
  MaxLinf <- MinLinf + 80
  
  init_f <- c(k = 0.10, zeta = 6, mu = log(180), sigma = 0.20)
  
  fit_f <- optim(
    par     = init_f,
    fn      = LL,
    dat     = dat_f,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.10, 1.01, log(120), 0.05),
    upper   = c(0.50, 200, log(260), 1.00),
    control = list(maxit = 300)
  )
  
  val_f <- fit_f$par
  se_f  <- get_se_k_Linf(fit_f, dat_f, MinLinf, MaxLinf)
  
  Linf_f <- exp(val_f[3] + 0.5 * val_f[4]^2)
  AIC_f  <- 2 * fit_f$value + 2 * length(val_f)
  
  # 5. MALE
  dat_m <- dat[dat$SEX == 2, ]
  MinLinf <- max(dat_m$PL + dat_m$INC) + 0.1
  MaxLinf <- MinLinf + 80
  
  init_m <- c(k = 0.010, zeta = 6, mu = log(180), sigma = 0.20)
  
  fit_m <- optim(
    par     = init_m,
    fn      = LL,
    dat     = dat_m,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.01, 1.01, log(120), 0.05),
    upper   = c(0.50, 200, log(260), 1.00),
    control = list(maxit = 300)
  )
  
  val_m <- fit_m$par
  se_m  <- get_se_k_Linf(fit_m, dat_m, MinLinf, MaxLinf)
  
  Linf_m <- exp(val_m[3] + 0.5 * val_m[4]^2)
  AIC_m  <- 2 * fit_m$value + 2 * length(val_m)
  
  # 6. FINAL TABLE (ONLY k & Linf SEs)
  table2_results <- data.frame(
    Sex      = c("Female", "Male"),
    k        = c(val_f[1], val_m[1]),
    k_SE     = c(se_f$se_k, se_m$se_k),
    Linf     = c(Linf_f, Linf_m),
    Linf_SE  = c(se_f$se_Linf, se_m$se_Linf),
    AIC      = c(AIC_f, AIC_m)
  )
  
  if (!dir.exists("results/tables"))
    dir.create("results/tables", recursive = TRUE)
  
  write.csv(
    table2_results,
    "results/tables/table2_bl_model.csv",
    row.names = FALSE
  )
  
  print(table2_results)
  invisible(table2_results)
}
test_table2_bl()
