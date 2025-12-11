test_table6joint <- function() {
  dat <- lobster
  
  ff <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li) {
    gLinf <- 0
    for (j in seq_along(Xi)) {
      a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
      b <- exp(-k * Ti[j]) * (zeta - 1)
      x <- Xi[j] / (Linf - Li[j])
      bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) -
        lgamma(a) - lgamma(b) -
        log(sigma) - 0.5 * log(2 * pi) -
        ((log(Linf) - miu)^2) / (2 * sigma^2) -
        log(Linf - Li[j])
      gLinf <- gLinf + bb
    }
    gLinf
  }
  
  gg <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
    res <- ff(Linf, k, zeta, miu, sigma, Xi, Ti, Li)
    Mi <- optimize(ff, c(MinLinf, MaxLinf),
                   k=k, zeta=zeta, miu=miu, sigma=sigma,
                   Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
    exp(res - Mi)
  }
  
  LL <- function(theta, dat, MinLinf, MaxLinf) {
    k <- theta[1]; zeta <- theta[2]; miu <- theta[3]; sigma <- theta[4]
    all.LL <- 0
    for (iid in unique(dat$LOBSTER)) {
      dati <- dat[dat$LOBSTER == iid, ]
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT/365.25
      Mi <- optimize(ff, c(MinLinf, MaxLinf),
                     k=k, zeta=zeta, miu=miu, sigma=sigma,
                     Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
      res.int <- integrate(gg, lower=MinLinf, upper=MaxLinf,
                           k=k, zeta=zeta, miu=miu, sigma=sigma,
                           Xi=Xi, Ti=Ti, Li=Li,
                           MinLinf=MinLinf, MaxLinf=MaxLinf)$value
      all.LL <- all.LL + log(res.int) + Mi
    }
    -all.LL
  }
  
  MinLinf <- 160; MaxLinf <- 200
  init_theta <- c(0.1, 4, 1, 2.88)
  num_params <- 4
  
  # Female
  dat_f <- subset(dat, SEX == 1)
  res_f <- optim(init_theta, LL, dat=dat_f,
                 MinLinf=MinLinf, MaxLinf=MaxLinf,
                 lower=c(0.01,0.1,0.1,0.1),
                 upper=c(1,100,10,10),
                 method="L-BFGS-B", control=list(maxit=1000))
  val_f <- res_f$par
  negLL_f <- res_f$value
  logLik_f <- -negLL_f
  AIC_f <- 2*negLL_f + 2*num_params
  Linf_f <- exp(val_f[3] + val_f[4]^2/2)
  
  # Male
  dat_m <- subset(dat, SEX == 2)
  res_m <- optim(init_theta, LL, dat=dat_m,
                 MinLinf=MinLinf, MaxLinf=MaxLinf,
                 lower=c(0.01,0.1,0.1,0.1),
                 upper=c(1,100,10,10),
                 method="L-BFGS-B", control=list(maxit=1000))
  val_m <- res_m$par
  negLL_m <- res_m$value
  logLik_m <- -negLL_m
  AIC_m <- 2*negLL_m + 2*num_params
  Linf_m <- exp(val_m[3] + val_m[4]^2/2)
  
  # Assemble results
  results_table <- data.frame(
    Sex       = c("Female","Male"),
    k         = c(val_f[1], val_m[1]),
    zeta      = c(val_f[2], val_m[2]),
    miu       = c(val_f[3], val_m[3]),
    sigma     = c(val_f[4], val_m[4]),
    E_Linf    = c(Linf_f, Linf_m),
    LogLik    = c(logLik_f, logLik_m),
    NegLogLik = c(negLL_f, negLL_m),
    AIC       = c(AIC_f, AIC_m),
    Converged = c(res_f$convergence==0, res_m$convergence==0)
  )
  
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive=TRUE)
  write.csv(results_table, "results/tables/table6_joint_models.csv", row.names=FALSE)

  invisible(results_table)
}
