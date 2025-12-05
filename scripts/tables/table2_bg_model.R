test_table2_bg <- function() {
  cat("Running test for table2 (BG model with Gamma prior)...\n")
  
  dat <- lobster
  
  # --- Helper functions ---
  ff_BG <- function(Linf, k, zeta, alpha, beta, Xi, Ti, Li) {
    if (Linf <= max(Li) + 1e-8) return(-1e12)
    log_prior <- (alpha - 1)*log(Linf) - (Linf/beta) -
      alpha*log(beta) - lgamma(alpha)
    out <- log_prior
    for (j in seq_along(Xi)) {
      denom <- Linf - Li[j]
      if (denom <= 1e-12) return(-1e12)
      x <- Xi[j] / denom
      if (x <= 1e-12 || x >= 1 - 1e-12) return(-1e12)
      a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
      b <- exp(-k * Ti[j]) * (zeta - 1)
      if (a <= 0 || b <= 0) return(-1e12)
      out <- out + dbeta(x, a, b, log = TRUE) - log(denom)
    }
    return(out)
  }
  
  kernel_BG <- function(Linf, Mi, k, zeta, alpha, beta, Xi, Ti, Li) {
    val <- ff_BG(Linf, k, zeta, alpha, beta, Xi, Ti, Li)
    exp(val - Mi)
  }
  
  safe_integrate <- function(fun, lower, upper, n=1200) {
    xs <- seq(lower, upper, length.out=n)
    vals <- sapply(xs, fun)
    dx <- (upper - lower)/(n-1)
    sum(vals) * dx
  }
  
  LL_BG <- function(theta, dat, MinLinf, MaxLinf) {
    k     <- theta[1]; zeta  <- theta[2]
    alpha <- theta[3]; beta  <- theta[4]
    if (k <= 0 || zeta <= 1.01 || alpha <= 0 || beta <= 0) return(1e12)
    LLtot <- 0
    for (id in unique(dat$LOBSTER)) {
      dati <- subset(dat, LOBSTER == id)
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      Mi <- optimize(
        function(Ls) ff_BG(Ls, k, zeta, alpha, beta, Xi, Ti, Li),
        interval = c(MinLinf, MaxLinf), maximum = TRUE
      )$objective
      val <- safe_integrate(
        function(Ls) kernel_BG(Ls, Mi, k, zeta, alpha, beta, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      if (!is.finite(val) || val <= 0) return(1e12)
      LLtot <- LLtot + log(val) + Mi
    }
    return(-LLtot)
  }
  
  # --- Female subset ---
  dat_female <- dat[dat$SEX == 1, ]
  MinLinf <- 180; MaxLinf <- 200
  init_theta_BG <- c(k=0.1, zeta=4, alpha=170, beta=1.0)
  
  res_female_BG <- optim(
    par     = init_theta_BG,
    fn      = LL_BG,
    dat     = dat_female,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.01, 1.02, 1e-4, 1e-4),
    upper   = c(2,    200,  500,  20),
    control = list(maxit = 1500)
  )
  
  val_female_BG <- res_female_BG$par
  AIC_female_BG <- 2 * res_female_BG$value + 2 * length(val_female_BG)
  Linf_female_BG <- val_female_BG[3] * val_female_BG[4]
  
  # --- Male subset ---
  dat_male <- dat[dat$SEX == 2, ]
  init_theta_BG_male <- c(k=0.1, zeta=4, alpha=10, beta=17)
  
  res_male_BG <- optim(
    par     = init_theta_BG_male,
    fn      = LL_BG,
    dat     = dat_male,
    MinLinf = MinLinf,
    MaxLinf = MaxLinf,
    method  = "L-BFGS-B",
    lower   = c(0.01, 0.1, 1, 1),
    upper   = c(1, 100, 500, 50),
    control = list(maxit = 5000)
  )
  
  val_male_BG <- res_male_BG$par
  logLik_male <- -res_male_BG$value
  AIC_male_BG <- 2 * length(val_male_BG) - 2 * logLik_male
  Linf_male_BG <- val_male_BG[3] * val_male_BG[4]
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex   = c("Female", "Male"),
    k     = c(val_female_BG[1], val_male_BG[1]),
    zeta  = c(val_female_BG[2], val_male_BG[2]),
    alpha = c(val_female_BG[3], val_male_BG[3]),
    beta  = c(val_female_BG[4], val_male_BG[4]),
    E_Linf = c(Linf_female_BG, Linf_male_BG),
    AIC    = c(AIC_female_BG, AIC_male_BG)
  )
  
  # --- Save to results/tables ---
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_bg_model.csv", row.names = FALSE)
    
  invisible(table2_results)
}
