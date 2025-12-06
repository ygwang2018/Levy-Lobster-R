test_table2_gg <- function() {
  
  dat <- lobster
  
  # --- Helper functions ---
  ff_GG <- function(Linf, k, alpha, a_L, b_L, Xi, Ti, Li) {
    if (!is.finite(Linf) || Linf <= max(Li) + 1e-10) return(-1e12)
    log_prior <- dgamma(Linf, shape = a_L, rate = b_L, log = TRUE)
    ll <- log_prior
    for (j in seq_along(Xi)) {
      mu <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
      if (mu <= 0 || Xi[j] <= 0) return(-1e12)
      shape <- alpha; scale <- mu / alpha
      if (shape <= 0 || scale <= 0) return(-1e12)
      ll <- ll + dgamma(Xi[j], shape = shape, scale = scale, log = TRUE)
    }
    return(ll)
  }
  
  kernel_GG <- function(Linf, Mi, k, alpha, a_L, b_L, Xi, Ti, Li) {
    val <- ff_GG(Linf, k, alpha, a_L, b_L, Xi, Ti, Li)
    exp(val - Mi)
  }
  
  safe_integrate <- function(fun, lower, upper, n = 100) {
    xs <- seq(lower, upper, length.out=n)
    vals <- sapply(xs, fun)
    dx <- (upper - lower)/(n-1)
    sum(vals) * dx
  }
  
  LL_GG <- function(theta, dat, MinLinf, MaxLinf) {
    k     <- theta[1]; alpha <- theta[2]
    a_L   <- theta[3]; b_L   <- theta[4]
    if (k <= 0 || alpha <= 0 || a_L <= 0 || b_L <= 0) return(1e12)
    LLtot <- 0
    for (id in unique(dat$LOBSTER)) {
      dati <- subset(dat, LOBSTER == id)
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      Mi <- optimize(
        function(Ls) ff_GG(Ls, k, alpha, a_L, b_L, Xi, Ti, Li),
        interval = c(MinLinf, MaxLinf), maximum = TRUE
      )$objective
      val <- safe_integrate(
        function(Ls) kernel_GG(Ls, Mi, k, alpha, a_L, b_L, Xi, Ti, Li),
        MinLinf, MaxLinf
      )
      if (!is.finite(val) || val <= 0) return(1e12)
      LLtot <- LLtot + log(val) + Mi
    }
    return(-LLtot)
  }
  
  # --- Female subset ---
  dat_f <- subset(dat, SEX == 1)
  MinLinf <- max(dat_f$PL) + 0.1; MaxLinf <- MinLinf + 80
  init_f <- c(k=0.15, alpha=10, a_L=180, b_L=1)
  res_f <- optim(init_f, LL_GG, dat=dat_f, MinLinf=MinLinf, MaxLinf=MaxLinf,
                 method="L-BFGS-B", lower=c(0.01,0.1,1e-3,1e-3),
                 upper=c(1.0,100,500,20), control=list(maxit=2000))
  val_f <- res_f$par; negLL_f <- res_f$value
  AIC_f <- 2*negLL_f + 2*length(val_f)
  Linf_f <- val_f[3]/val_f[4]
  
  # --- Male subset ---
  dat_m <- subset(dat, SEX == 2)
  MinLinf <- max(dat_m$PL) + 0.1; MaxLinf <- MinLinf + 80
  init_m <- c(k=0.15, alpha=10, a_L=160, b_L=1)
  res_m <- optim(init_m, LL_GG, dat=dat_m, MinLinf=MinLinf, MaxLinf=MaxLinf,
                 method="L-BFGS-B", lower=c(0.01,0.1,1e-3,1e-3),
                 upper=c(1.0,100,500,20), control=list(maxit=2000))
  val_m <- res_m$par; negLL_m <- res_m$value
  AIC_m <- 2*negLL_m + 2*length(val_m)
  Linf_m <- val_m[3]/val_m[4]
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex    = c("Female", "Male"),
    k      = c(val_f[1], val_m[1]),
    alpha  = c(val_f[2], val_m[2]),
    a_L    = c(val_f[3], val_m[3]),
    b_L    = c(val_f[4], val_m[4]),
    E_Linf = c(Linf_f, Linf_m),
    AIC    = c(AIC_f, AIC_m)
  )
  
  # --- Save to results/tables ---
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_gg_model.csv", row.names = FALSE)
    
  invisible(table2_results)
}
