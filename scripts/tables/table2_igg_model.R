  test_table2_igg <- function() {
  
  dat <- lobster
  
  # --- Helper functions ---
  log_MI_gamma <- function(Delta, Linf, lambda, k, Li, Ti) {
    Iij <- (Linf - Li) * (1 - exp(-k * Ti))
    shape <- lambda * Iij
    rate  <- lambda
    if (Iij <= 0 || Delta <= 0 || shape <= 0 || rate <= 0) return(-1e12)
    shape*log(rate) + (shape-1)*log(Delta) - rate*Delta - lgamma(shape)
  }
  
  log_prior_Linf <- function(Linf, alpha, beta) {
    dgamma(Linf, shape=alpha, rate=beta, log=TRUE)
  }
  
  loglik_per_individual <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta) {
    s <- log_prior_Linf(Linf, alpha, beta)
    for (j in seq_along(Xi)) {
      s <- s + log_MI_gamma(Xi[j], Linf, lambda, k, Li[j], Ti[j])
    }
    s
  }
  
  integrand <- function(Linf, Mi, Xi, Li, Ti, lambda, k, alpha, beta) {
    val <- loglik_per_individual(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
    exp(val - Mi)   # stabilised
  }
  
  LL <- function(theta, dat, MinLinf, MaxLinf) {
    lambda <- theta[1]; k <- theta[2]; alpha <- theta[3]; beta <- theta[4]
    if (lambda <= 0 || k <= 0 || alpha <= 0 || beta <= 0) return(1e12)
    LL_total <- 0
    for (iid in unique(dat$LOBSTER)) {
      dati <- dat[dat$LOBSTER == iid, ]
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      if (any(Li + Xi >= MaxLinf)) next
      Mi <- optimize(function(Ls) loglik_per_individual(Ls, Xi, Li, Ti, lambda, k, alpha, beta),
                     interval=c(MinLinf, MaxLinf), maximum=TRUE)$objective
      val <- try(integrate(integrand, lower=MinLinf, upper=MaxLinf,
                           Mi=Mi, Xi=Xi, Li=Li, Ti=Ti,
                           lambda=lambda, k=k, alpha=alpha, beta=beta,
                           rel.tol=1e-4)$value, silent=TRUE)
      if (inherits(val,"try-error") || val <= 0 || !is.finite(val)) return(1e12)
      LL_total <- LL_total + log(val) + Mi
    }
    -LL_total
  }
  
  # --- Female subset ---
  dat_f <- subset(dat, SEX == 1)
  MinLinf_f <- max(dat_f$PL) + 0.1; MaxLinf_f <- MinLinf_f + 80
  start_f <- c(lambda=0.1, k=0.2, alpha=5, beta=0.05)
  fit_f <- optim(start_f, LL, dat=dat_f, MinLinf=MinLinf_f, MaxLinf=MaxLinf_f,
                 method="L-BFGS-B", lower=c(0.01,0.01,0.10,0.01), upper=c(1,0.5,10,100),
                 control=list(maxit=2000))
  par_f <- fit_f$par; negLL_f <- fit_f$value
  AIC_f <- 2*negLL_f + 2*length(par_f)
  E_Linf_f <- par_f["alpha"]/par_f["beta"]
  
  # --- Male subset ---
  dat_m <- subset(dat, SEX == 2)
  MinLinf_m <- max(dat_m$PL) + 0.1; MaxLinf_m <- MinLinf_m + 80
  start_m <- c(lambda=0.1, k=0.2, alpha=5, beta=0.05)
  fit_m <- optim(start_m, LL, dat=dat_m, MinLinf=MinLinf_m, MaxLinf=MaxLinf_m,
                 method="L-BFGS-B", lower=c(0.01,0.01,0.10,0.01), upper=c(1,0.5,10,100),
                 control=list(maxit=2000))
  par_m <- fit_m$par; negLL_m <- fit_m$value
  AIC_m <- 2*negLL_m + 2*length(par_m)
  E_Linf_m <- par_m["alpha"]/par_m["beta"]
  
  # --- Combine and save ---
  table_results <- data.frame(
    Sex    = c("Female","Male"),
    lambda = c(par_f["lambda"], par_m["lambda"]),
    k      = c(par_f["k"],      par_m["k"]),
    alpha  = c(par_f["alpha"],  par_m["alpha"]),
    beta   = c(par_f["beta"],   par_m["beta"]),
    E_Linf = c(E_Linf_f,        E_Linf_m),
    AIC    = c(AIC_f,           AIC_m)
  )
  
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table_results, "results/tables/table2_igg_model.csv", row.names = FALSE)
    
  invisible(table_results)
}
