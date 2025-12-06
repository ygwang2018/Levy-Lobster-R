test_table2_gl <- function() {
  
  dat <- lobster
  
  # --- Helper functions ---
  ff_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li) {
    log_prior <- dlnorm(Linf, meanlog = miu, sdlog = sigma, log = TRUE)
    gLinf <- log_prior
    for (j in seq_along(Xi)) {
      mu    <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
      shape <- alpha
      scale <- mu / alpha
      gLinf <- gLinf + dgamma(Xi[j], shape = shape, scale = scale, log = TRUE)
    }
    return(gLinf)
  }
  
  gg_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
    res <- ff_GL(Linf, k, alpha, miu, sigma, Xi, Ti, Li)
    Mi  <- optimize(ff_GL, c(MinLinf, MaxLinf),
                    k=k, alpha=alpha, miu=miu, sigma=sigma,
                    Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
    exp(res - Mi)
  }
  
  LL_GL <- function(theta, dat, MinLinf, MaxLinf) {
    k <- theta[1]; alpha <- theta[2]; miu <- theta[3]; sigma <- theta[4]
    all.LL <- 0
    for (iid in unique(dat$LOBSTER)) {
      dati <- dat[dat$LOBSTER == iid, ]
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      if (any(Li + Xi >= MaxLinf)) next
      Mi <- optimize(ff_GL, c(MinLinf, MaxLinf),
                     k=k, alpha=alpha, miu=miu, sigma=sigma,
                     Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
      res.int <- tryCatch({
        integrate(gg_GL, lower=MinLinf, upper=MaxLinf,
                  k=k, alpha=alpha, miu=miu, sigma=sigma,
                  Xi=Xi, Ti=Ti, Li=Li,
                  MinLinf=MinLinf, MaxLinf=MaxLinf,
                  rel.tol=1e-4)$value
      }, error=function(e) NA)
      if (is.na(res.int) || res.int <= 0) next
      all.LL <- all.LL + log(res.int) + Mi
    }
    return(-all.LL)
  }
  
  MinLinf <- 160; MaxLinf <- 200
  init_theta <- c(0.1, 4, 5, 0.5)
  num_params <- 4
  
  # --- Female subset ---
  dat_f <- subset(dat, SEX == 1)
  res_f <- optim(init_theta, LL_GL, dat=dat_f,
                 MinLinf=MinLinf, MaxLinf=MaxLinf,
                 lower=c(0.01,0.1,0.1,0.01),
                 upper=c(1,100,10,10),
                 method="L-BFGS-B", control=list(maxit=1000))
  val_f <- res_f$par; negLL_f <- res_f$value
  AIC_f <- 2*num_params + 2*negLL_f
  Linf_f <- exp(val_f[3] + (val_f[4]^2)/2)
  
  # --- Male subset ---
  dat_m <- subset(dat, SEX == 2)
  res_m <- optim(init_theta, LL_GL, dat=dat_m,
                 MinLinf=MinLinf, MaxLinf=MaxLinf,
                 lower=c(0.01,0.1,0.1,0.01),
                 upper=c(1,100,10,10),
                 method="L-BFGS-B", control=list(maxit=1000))
  val_m <- res_m$par; negLL_m <- res_m$value
  AIC_m <- 2*num_params + 2*negLL_m
  Linf_m <- exp(val_m[3] + (val_m[4]^2)/2)
  
  # --- Combine results ---
  table2_results <- data.frame(
    Sex    = c("Female", "Male"),
    k      = c(val_f[1], val_m[1]),
    alpha  = c(val_f[2], val_m[2]),
    miu    = c(val_f[3], val_m[3]),
    sigma  = c(val_f[4], val_m[4]),
    E_Linf = c(Linf_f, Linf_m),
    AIC    = c(AIC_f, AIC_m)
  )
  
  # --- Save to results/tables ---
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table2_results, "results/tables/table2_gl_model.csv", row.names = FALSE)
  
  
  invisible(table2_results)
}
