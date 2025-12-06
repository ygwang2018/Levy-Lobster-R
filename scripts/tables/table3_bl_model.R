test_table3_bl <- function() {
  library(dplyr)

  max_over_Linf_safe <- function(fun, lower, upper) {
    opt <- try(optimize(fun, interval = c(lower, upper), maximum = TRUE), silent = TRUE)
    obj1 <- if (inherits(opt, "try-error")) -1e10 else opt$objective
    vals <- c(obj1, fun(lower), fun(upper), fun((lower + upper) / 2))
    max(vals[is.finite(vals)])
  }

  safe_integrate <- function(fun, lower, upper, n = 300) {
    xs <- seq(lower, upper, length.out = n)
    vals <- sapply(xs, function(x) {
      v <- fun(x)
      if (!is.finite(v)) 0 else v
    })
    dx <- (upper - lower) / (n - 1)
    sum(vals) * dx
  }

  ff_ind <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
    if (!is.finite(Linf) || Linf <= max(Li) + 1e-8) return(-1e10)
    logPrior <- -log(Linf) - log(sigma) - 0.5*log(2*pi) - ((log(Linf) - mu)^2 / (2*sigma^2))
    out <- logPrior
    for (j in seq_along(Xi)) {
      denom <- Linf - Li[j]; if (denom <= 1e-8) return(-1e10)
      x <- Xi[j] / denom;   if (x <= 1e-8 || x >= 1 - 1e-8) return(-1e10)
      a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
      b <- exp(-k * Ti[j]) * (zeta - 1)
      if (a <= 1e-8 || b <= 1e-8) return(-1e-10)
      out <- out + dbeta(x, a, b, log = TRUE) - log(denom)
    }
    out
  }

  kernel_norm <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
    val <- ff_ind(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
    exp(val - Mi)
  }

  LL <- function(theta, dat, MinLinf, MaxLinf) {
    k <- theta[1]; zeta <- theta[2]; mu <- theta[3]; sigma <- theta[4]
    if (k <= 0 || zeta <= 1.01 || sigma <= 0) return(1e12)
    LLtot <- 0
    for (id in unique(dat$LOBSTER)) {
      dati <- subset(dat, LOBSTER == id)
      Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
      Mi <- max_over_Linf_safe(function(Ls) ff_ind(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
                               MinLinf, MaxLinf)
      val <- safe_integrate(function(Ls) kernel_norm(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
                            MinLinf, MaxLinf)
      if (val <= 0 || !is.finite(val)) return(1e12)
      LLtot <- LLtot + log(val) + Mi
    }
    -LLtot
  }

  # Female
  dat_f <- lobster %>% filter(SEX == 1)
  MinLinf_f <- max(dat_f$PL + dat_f$INC) + 0.1
  MaxLinf_f <- MinLinf_f + 80
  init_f <- c(k = 0.10, zeta = 6, mu = log(180), sigma = 0.20)
  fit_f <- optim(par = init_f, fn = LL, dat = dat_f, MinLinf = MinLinf_f, MaxLinf = MaxLinf_f,
                 method = "L-BFGS-B",
                 lower = c(0.1, 1.01, log(120), 0.05),
                 upper = c(0.5, 200, log(260), 1.0),
                 control = list(maxit = 300))
  par_f <- fit_f$par
  AIC_f <- 2 * fit_f$value + 2 * length(par_f)
  E_Linf_f <- exp(par_f[3] + par_f[4]^2 / 2)

  # Male
  dat_m <- lobster %>% filter(SEX == 2)
  MinLinf_m <- max(dat_m$PL + dat_m$INC) + 0.1
  MaxLinf_m <- MinLinf_m + 80
  init_m <- c(k = 0.010, zeta = 6, mu = log(180), sigma = 0.20)
  fit_m <- optim(par = init_m, fn = LL, dat = dat_m, MinLinf = MinLinf_m, MaxLinf = MaxLinf_m,
                 method = "L-BFGS-B",
                 lower = c(0.01, 1.01, log(120), 0.05),
                 upper = c(0.5, 200, log(260), 1.0),
                 control = list(maxit = 300))
  par_m <- fit_m$par
  AIC_m <- 2 * fit_m$value + 2 * length(par_m)
  E_Linf_m <- exp(par_m[3] + par_m[4]^2 / 2)

  table_bl <- data.frame(
    Sex    = c("Female", "Male"),
    k      = c(par_f[1], par_m[1]),
    zeta   = c(par_f[2], par_m[2]),
    mu     = c(par_f[3], par_m[3]),
    sigma  = c(par_f[4], par_m[4]),
    AIC    = c(AIC_f, AIC_m),
    E_Linf = c(E_Linf_f, E_Linf_m)
  )

  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(table_bl, "results/tables/table3_bl_model.csv", row.names = FALSE)
 
  invisible(table_bl)
}
