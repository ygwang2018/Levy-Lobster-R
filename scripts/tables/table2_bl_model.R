test_table2_bl <- function() {
  cat("Running test for table2 (BL model)...\n")
  
  dat <- lobster
  
  # --- Helper functions ---
  ff_ind <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
    if (!is.finite(Linf) || Linf <= max(Li) + 1e-8) return(-1e10)
    logPrior <- -log(Linf) - log(sigma) - 0.5*log(2*pi) -
      ((log(Linf) - mu)^2 / (2*sigma^2))
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
    return(out)
  }
  
  kernel_norm <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
    val <- ff_ind(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
    exp(val - Mi)
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
    dx <- (upper - lower)/(n-1)
    sum(vals) * dx
  }
  
  LL <- function(theta, dat, MinLinf, MaxLinf) {
    k     <- theta[1]; zeta  <- theta[2]
    mu    <- theta[3]; sigma <- theta[4]
    if (k <= 0 || zeta <= 1.01 || sigma <= 0) return(1e12)
    LLtot <- 0
    for (id in unique(dat$LOBSTER)) {
      dati <- subset(dat, LOBSTER == id)
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      Mi <- max_over_Linf(
        function(Ls) ff_ind(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      val <- safe_integrate(
        function(Ls) kernel_norm(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      if (val <= 0 || !is.finite(val)) return(1e12)
      LLtot <- LLtot + log(val) + Mi
    }
    return(-LLtot)
  }
  
  # --- Female subset ---
  dat_f <- subset(dat, SEX == 1)
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
    lower   = c(0.1, 1.01, log(120), 0.05),
    upper   = c(0.5, 200, log(260), 1.0),
    control = list(maxit = 300)
  )
  
  val_f <- fit_f$par
  negLL_f <- fit_f$value
  AIC_f <- 2*negLL_f + 2*length(val_f)
  sigma_f <- val_f[4]
  Linf_f <- exp(val_f[3] + sigma_f^2/2)
  
  # --- Male subset ---
  dat_m <- subset(dat, SEX == 2)
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
    upper   = c(0.5, 200, log(260), 1.0),
    control = list(maxit = 300)
  )
  
  val_m <- fit_m$par
  negLL_m <- fit_m$value
  AIC_m <- 2*negLL_m + 2*length(val_m)
  sigma_m <- val_m[4]
  Linf_m <- exp(val_m[3] + sigma_m^2/2)
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex   = c("Female", "Male"),
    k     = c(val_f[1], val_m[1]),
    zeta  = c(val_f[2], val_m[2]),
    mu    = c(val_f[3], val_m[3]),
    sigma = c(val_f[4], val_m[4]),
    E_Linf = c(Linf_f, Linf_m),
    AIC    = c(AIC_f, AIC_m)
  )
  
  # --- Save to results/tables ---
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_bl_model.csv", row.names = FALSE)
  
  invisible(table2_results)
}
