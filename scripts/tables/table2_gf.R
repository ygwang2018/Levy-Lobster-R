test_table2_gf <- function() {
  cat("Running test for table2 (Gamma–Linf fixed model)...\n")
  
  dat <- lobster
  
  LL1 <- function(d, L, I, T) {
    Linf   <- d[1]
    lambda <- d[2]
    k      <- d[3]
    
    # Constraints
    if (Linf <= max(L) + 1e-6) return(1e10)
    if (lambda <= 0 || k <= 0) return(1e10)
    
    mu <- (Linf - L) * (1 - exp(-k * T))
    if (any(mu <= 0) || any(I <= 0)) return(1e10)
    
    aa <- lambda * mu
    if (any(aa <= 0)) return(1e10)
    
    ll <- aa * log(lambda) +
      (aa - 1) * log(I) -
      lambda * I -
      lgamma(aa)
    
    return(-sum(ll))
  }
  
  # --- Female subset ---
  f <- dat[dat$SEX == 1, ]
  T <- f$INT / 365.25; I <- f$INC; L <- f$PL
  res_f <- optim(
    par    = c(180, 0.5, 0.01),
    fn     = LL1,
    L      = L, I = I, T = T,
    method = "L-BFGS-B",
    lower  = c(160, 0.1, 0.01),
    upper  = c(300, 1, 0.5)
  )
  val_f <- res_f$par
  AIC_f <- 2*res_f$value + 2*length(val_f)
  
  # --- Male subset ---
  f <- dat[dat$SEX == 2, ]
  T <- f$INT / 365.25; I <- f$INC; L <- f$PL
  res_m <- optim(
    par    = c(180, 0.5, 0.01),
    fn     = LL1,
    L      = L, I = I, T = T,
    method = "L-BFGS-B",
    lower  = c(180, 0.1, 0.01),
    upper  = c(300, 1, 0.5)
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
  write.csv(table2_results, "results/tables/table2_gf.csv", row.names = FALSE)
    
  invisible(table2_results)
}
