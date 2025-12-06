test_table2_ig <- function() {
  
  dat <- lobster
  
  # --- Helper functions ---
  log_dinvgauss <- function(x, mu, lambda) {
    0.5 * (log(lambda) - log(2*pi) - 3*log(x)) -
      lambda * (x - mu)^2 / (2 * mu^2 * x)
  }
  
  LL_IG <- function(d, L, I, T) {
    Linf   <- d[1]; lambda <- d[2]; k <- d[3]
    if (Linf <= max(L) + 1e-6 || lambda <= 0 || k <= 0) return(1e10)
    mu <- (Linf - L) * (1 - exp(-k * T))
    if (any(mu <= 0) || any(I <= 0)) return(1e10)
    ll <- log_dinvgauss(I, mu, lambda)
    return(-sum(ll))
  }
  
  # --- Female subset ---
  f <- dat[dat$SEX == 1, ]
  T <- f$INT / 365.25; I <- f$INC; L <- f$PL
  res_f <- optim(
    par    = c(180, 1, 0.1),
    fn     = LL_IG,
    L      = L, I = I, T = T,
    method = "L-BFGS-B",
    lower  = c(100, 0.01, 0.01),
    upper  = c(300, 10, 0.4),
    control = list(maxit = 2000)
  )
  val_f <- res_f$par
  AIC_f <- 2*res_f$value + 2*length(val_f)
  
  # --- Male subset ---
  f <- dat[dat$SEX == 2, ]
  T <- f$INT / 365.25; I <- f$INC; L <- f$PL
  res_m <- optim(
    par    = c(160, 1, 0.1),
    fn     = LL_IG,
    L      = L, I = I, T = T,
    method = "L-BFGS-B",
    lower  = c(160, 0.01, 0.01),
    upper  = c(180, 10, 1),
    control = list(maxit = 2000)
  )
  val_m <- res_m$par
  AIC_m <- 2*res_m$value + 2*length(val_m)
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex    = c("Female", "Male"),
    Linf   = c(val_f[1], val_m[1]),
    lambda = c(val_f[2], val_m[2]),
    k      = c(val_f[3], val_m[3]),
    AIC    = c(AIC_f, AIC_m)
  )
  
  # --- Save to results/tables ---
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_ig_model.csv", row.names = FALSE)
    
  invisible(table2_results)
}
