library(dplyr)

test_table2_bf <- function() {
  
  dat <- lobster
  
  ## --- Female subset ---
  f <- dat[dat$SEX == 1, ]
  T <- f$INT / 365.25
  I <- f$INC
  L <- f$PL
  N <- nrow(f)
  
  LL1 <- function(d) {
    Linf <- d[1]; k <- d[2]; zeta <- d[3]
    if (Linf <= max(L) + 1e-6 || k <= 0 || zeta <= 1) return(1e10)
    BG <- 0
    for (i in 1:N) {
      if (I[i] >= Linf - L[i]) return(1e10)
      alpha <- (1 - exp(-k*T[i])) * (zeta-1)
      beta  <- exp(-k*T[i]) * (zeta-1)
      if (alpha <= 0 || beta <= 0) return(1e10)
      x <- I[i] / (Linf - L[i])
      BG <- BG +
        lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) +
        (alpha-1)*log(x) + (beta-1)*log(1-x) -
        log(Linf - L[i])
    }
    -BG
  }
  
  res_f <- optim(c(200, 0.1, 50), LL1)$par
  Linf_f <- res_f[1]; k_f <- res_f[2]; zeta_f <- res_f[3]
  
  ## --- Male subset ---
  m <- dat[dat$SEX == 2, ]
  T <- m$INT / 365.25
  I <- m$INC
  L <- m$PL
  N <- nrow(m)
  
  res_m <- optim(c(200, 0.1, 50), LL1)$par
  Linf_m <- res_m[1]; k_m <- res_m[2]; zeta_m <- res_m[3]
  
  ## --- Combine results ---
  table2_results <- data.frame(
    Sex   = c("Female", "Male"),
    Linf  = c(Linf_f, Linf_m),
    k     = c(k_f, k_m),
    zeta  = c(zeta_f, zeta_m)
  )
  
  # Save to results/tables
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_bf.csv", row.names = FALSE)
    
  invisible(table2_results)
}

