test_table2_gg <- function() {
  
  dat <- lobster
  
  # 1. Log-posterior for GG model
  ff_GG <- function(Linf, k, alpha, a_L, b_L, Xi, Ti, Li) {
    
    if (!is.finite(Linf) || Linf <= max(Li) + 1e-10) return(-1e12)
    
    # Gamma prior on Linf
    ll <- dgamma(Linf, shape = a_L, rate = b_L, log = TRUE)
    
    for (j in seq_along(Xi)) {
      
      mu <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
      if (mu <= 0 || Xi[j] <= 0) return(-1e12)
      
      shape <- alpha
      scale <- mu / alpha
      
      if (shape <= 0 || scale <= 0) return(-1e12)
      
      ll <- ll + dgamma(Xi[j], shape = shape, scale = scale, log = TRUE)
    }
    
    ll
  }
  
  kernel_GG <- function(Linf, Mi, k, alpha, a_L, b_L, Xi, Ti, Li) {
    exp(ff_GG(Linf, k, alpha, a_L, b_L, Xi, Ti, Li) - Mi)
  }
  
  safe_integrate <- function(fun, lower, upper, n = 100) {
    xs <- seq(lower, upper, length.out = n)
    vals <- sapply(xs, fun)
    dx <- (upper - lower) / (n - 1)
    sum(vals) * dx
  }
  
  # 2. Marginal log-likelihood
  LL_GG <- function(theta, dat, MinLinf, MaxLinf) {
    
    k     <- theta[1]
    alpha <- theta[2]
    a_L   <- theta[3]
    b_L   <- theta[4]
    
    if (k <= 0 || alpha <= 0 || a_L <= 0 || b_L <= 0)
      return(1e12)
    
    LLtot <- 0
    
    for (id in unique(dat$LOBSTER)) {
      
      dati <- dat[dat$LOBSTER == id, ]
      Xi <- dati$INC
      Li <- dati$PL
      Ti <- dati$INT / 365.25
      
      Mi <- optimize(
        function(Ls) ff_GG(Ls, k, alpha, a_L, b_L, Xi, Ti, Li),
        interval = c(MinLinf, MaxLinf),
        maximum = TRUE
      )$objective
      
      val <- safe_integrate(
        function(Ls) kernel_GG(Ls, Mi, k, alpha, a_L, b_L, Xi, Ti, Li),
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
      opt$par, fn = LL_GG,
      dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf
    )
    
    vcov <- tryCatch(
      solve(hess),
      error = function(e) matrix(NA, 4, 4)
    )
    
    # SE for k
    se_k <- sqrt(vcov[1, 1])
    
    # Delta method for E[Linf] = a_L / b_L
    a_L <- opt$par[3]
    b_L <- opt$par[4]
    
    grad <- c(
      d_a =  1 / b_L,
      d_b = -a_L / b_L^2
    )
    
    var_Linf <- grad %*% vcov[3:4, 3:4] %*% grad
    se_Linf  <- sqrt(var_Linf)
    
    list(se_k = se_k, se_Linf = se_Linf)
  }
  
  # 4. FEMALE
  dat_f <- dat[dat$SEX == 1, ]
  MinLinf <- max(dat_f$PL) + 0.1
  MaxLinf <- MinLinf + 80
  
  init_f <- c(k = 0.15, alpha = 10, a_L = 180, b_L = 1)
  
  res_f <- optim(
    par     = init_f,
    fn      = LL_GG,
    dat     = dat_f,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.01, 0.1, 1e-3, 1e-3),
    upper   = c(1.0, 100, 500, 20),
    control = list(maxit = 2000)
  )
  
  val_f <- res_f$par
  Linf_f <- val_f[3] / val_f[4]
  se_f   <- get_se_k_Linf(res_f, dat_f, MinLinf, MaxLinf)
  
  AIC_f <- 2 * res_f$value + 2 * length(val_f)
  
  # 5. MALE
  dat_m <- dat[dat$SEX == 2, ]
  MinLinf <- max(dat_m$PL) + 0.1
  MaxLinf <- MinLinf + 80
  
  init_m <- c(k = 0.15, alpha = 10, a_L = 160, b_L = 1)
  
  res_m <- optim(
    par     = init_m,
    fn      = LL_GG,
    dat     = dat_m,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.01, 0.1, 1e-3, 1e-3),
    upper   = c(1.0, 100, 500, 20),
    control = list(maxit = 2000)
  )
  
  val_m <- res_m$par
  Linf_m <- val_m[3] / val_m[4]
  se_m   <- get_se_k_Linf(res_m, dat_m, MinLinf, MaxLinf)
  
  AIC_m <- 2 * res_m$value + 2 * length(val_m)
  
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
    "results/tables/table2_gg_model.csv",
    row.names = FALSE
  )
  
  print(table2_results)
  invisible(table2_results)
}
test_table2_gg()
