test_table3_bg <- function() {
  
  dat <- lobster
  
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
    out
  }
  
  kernel_BG <- function(Linf, Mi, k, zeta, alpha, beta, Xi, Ti, Li) {
    val <- ff_BG(Linf, k, zeta, alpha, beta, Xi, Ti, Li)
    exp(val - Mi)
  }
  
  safe_integrate <- function(fun, lower, upper, n=100) {
    xs <- seq(lower, upper, length.out=n)
    vals <- sapply(xs, function(x) ifelse(is.finite(fun(x)), fun(x), 0))
    dx <- (upper - lower)/(n-1)
    sum(vals) * dx
  }
  
  LL_BG <- function(theta, dat, MinLinf, MaxLinf) {
    k     <- theta[1]; zeta <- theta[2]; alpha <- theta[3]; beta <- theta[4]
    if (k <= 0 || zeta <= 1.01 || alpha <= 0 || beta <= 0) return(1e12)
    LLtot <- 0
    for (id in unique(dat$LOBSTER)) {
      dati <- subset(dat, LOBSTER == id)
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      Mi <- optimize(function(Ls) ff_BG(Ls, k, zeta, alpha, beta, Xi, Ti, Li),
                     interval=c(MinLinf, MaxLinf), maximum=TRUE)$objective
      val <- safe_integrate(function(Ls) kernel_BG(Ls, Mi, k, zeta, alpha, beta, Xi, Ti, Li),
                            MinLinf, MaxLinf)
      if (!is.finite(val) || val <= 0) return(1e12)
      LLtot <- LLtot + log(val) + Mi
    }
    -LLtot
  }
  
  # --- Female subset ---
  dat_f <- subset(dat, SEX == 1)
  MinLinf <- 180; MaxLinf <- 200
  init_f <- c(k=0.1, zeta=4, alpha=170, beta=1.0)
  res_f <- optim(init_f, LL_BG, dat=dat_f, MinLinf=MinLinf, MaxLinf=MaxLinf,
                 method="L-BFGS-B", lower=c(0.01,1.02,1e-4,1e-4),
                 upper=c(2,200,500,20), control=list(maxit=1500))
  val_f <- res_f$par; AIC_f <- 2*res_f$value + 2*length(val_f)
  Linf_f <- val_f[3]*val_f[4]
  
  # --- Male subset ---
  dat_m <- subset(dat, SEX == 2)
  init_m <- c(0.1,4,10,17)
  res_m <- optim(init_m, LL_BG, dat=dat_m, MinLinf=MinLinf, MaxLinf=MaxLinf,
                 method="L-BFGS-B", lower=c(0.01,0.1,1,1),
                 upper=c(1,100,500,50), control=list(maxit=5000))
  val_m <- res_m$par; AIC_m <- 2*res_m$value + 2*length(val_m)
  Linf_m <- val_m[3]*val_m[4]
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex    = c("Female","Male"),
    k      = c(val_f[1], val_m[1]),
    zeta   = c(val_f[2], val_m[2]),
    alpha  = c(val_f[3], val_m[3]),
    beta   = c(val_f[4], val_m[4]),
    E_Linf = c(Linf_f, Linf_m),
    AIC    = c(AIC_f, AIC_m)
  )
  
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_bg_model.csv", row.names = FALSE)
  
  invisible(table2_results)
}
