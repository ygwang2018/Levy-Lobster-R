library(dplyr)

############################################################
#  GENERIC SAFE FUNCTIONS
############################################################

max_over_Linf_safe <- function(fun, lower, upper) {
  opt <- try(optimize(fun, interval = c(lower, upper), maximum = TRUE), silent = TRUE)
  if (inherits(opt, "try-error")) {
    obj1 <- -1e10
  } else {
    obj1 <- opt$objective
  }
  vals <- c(obj1, fun(lower), fun(upper), fun((lower + upper) / 2))
  max(vals[is.finite(vals)])
}

safe_integrate <- function(fun, lower, upper, n = 120) {
  xs <- seq(lower, upper, length.out = n)
  vals <- sapply(xs, fun)
  vals[!is.finite(vals)] <- 0
  dx <- (upper - lower) / (n - 1)
  sum(vals) * dx
}

############################################################
#  FEMALE MODEL
############################################################

dat_f <- lobster %>% filter(SEX == 1)
MinLinf_f <- max(dat_f$PL + dat_f$INC) + 0.1
MaxLinf_f <- MinLinf_f + 80

ff_ind_f <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {

  if (!is.finite(Linf) || Linf <= max(Li) + 1e-8)
    return(-1e10)

  logPrior <- -log(Linf) - log(sigma) - 0.5*log(2*pi) -
    ((log(Linf) - mu)^2 / (2*sigma^2))

  out <- logPrior

  for (j in seq_along(Xi)) {
    denom <- Linf - Li[j]
    if (denom <= 1e-8) return(-1e10)

    x <- Xi[j] / denom
    if (x <= 1e-8 || x >= 1 - 1e-8) return(-1e10)

    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    if (a <= 1e-8 || b <= 1e-8) return(-1e10)

    out <- out + dbeta(x, a, b, log = TRUE) - log(denom)
  }

  out
}

kernel_norm_f <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
  val <- ff_ind_f(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
  exp(val - Mi)
}

LL_f <- function(theta, dat, MinLinf, MaxLinf) {

  k     <- theta[1]
  zeta  <- theta[2]
  mu    <- theta[3]
  sigma <- theta[4]

  if (k <= 0 || zeta <= 1.01 || sigma <= 0)
    return(1e12)

  LLtot <- 0

  for (id in unique(dat$LOBSTER)) {

    dati <- subset(dat, LOBSTER == id)

    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25

    Mi <- max_over_Linf_safe(
      function(Ls) ff_ind_f(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )

    val <- safe_integrate(
      function(Ls) kernel_norm_f(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )

    if (val <= 0 || !is.finite(val))
      return(1e12)

    LLtot <- LLtot + log(val) + Mi
  }

  -LLtot
}

# fit female
init_f <- c(k = 0.10, zeta = 6, mu = log(180), sigma = 0.20)

fit_f <- optim(
  par     = init_f,
  fn      = LL_f,
  dat     = dat_f,
  MinLinf = MinLinf_f,
  MaxLinf = MaxLinf_f,
  method  = "L-BFGS-B",
  lower   = c(0.1, 1.01, log(120), 0.05),
  upper   = c(0.5, 200, log(260), 1.0),
  control = list(maxit = 300)
)

negLL_f <- fit_f$value
p <- length(init_f)
AIC_f <- 2 * negLL_f + 2 * p

# CI for female
par_f <- fit_f$par
H_f <- optimHess(par_f, LL_f, dat = dat_f, MinLinf = MinLinf_f, MaxLinf = MaxLinf_f)
cov_f <- solve(H_f)
se_f  <- sqrt(diag(cov_f))

# CI for k_f
CI_k_f <- par_f[1] + c(-1.96, 1.96) * se_f[1]

# CI for E[Linf]
Linf_f <- exp(par_f[3] + par_f[4]^2 / 2)
grad_Linf_f <- c(0,0, Linf_f, par_f[4] * Linf_f)
var_Linf_f <- t(grad_Linf_f) %*% cov_f %*% grad_Linf_f
se_Linf_f <- sqrt(var_Linf_f)
CI_Linf_f <- Linf_f + c(-1.96, 1.96) * se_Linf_f


############################################################
#  MALE MODEL (same structure)
############################################################

dat_m <- lobster %>% filter(SEX == 2)
MinLinf_m <- max(dat_m$PL + dat_m$INC) + 0.1
MaxLinf_m <- MinLinf_m + 80

# identical functions but renamed to *_m
ff_ind_m <- ff_ind_f
kernel_norm_m <- kernel_norm_f
LL_m <- LL_f

init_m <- c(k = 0.010, zeta = 6, mu = log(180), sigma = 0.20)

fit_m <- optim(
  par     = init_m,
  fn      = LL_m,
  dat     = dat_m,
  MinLinf = MinLinf_m,
  MaxLinf = MaxLinf_m,
  method  = "L-BFGS-B",
  lower   = c(0.01, 1.01, log(120), 0.05),
  upper   = c(0.5, 200, log(260), 1.0),
  control = list(maxit = 300)
)

negLL_m <- fit_m$value
AIC_m <- 2 * negLL_m + 2 * length(init_m)

# CI male
par_m <- fit_m$par
H_m <- optimHess(par_m, LL_m, dat = dat_m, MinLinf = MinLinf_m, MaxLinf = MaxLinf_m)
cov_m <- solve(H_m)
se_m  <- sqrt(diag(cov_m))

CI_k_m <- par_m[1] + c(-1.96,1.96) * se_m[1]

Linf_m <- exp(par_m[3] + par_m[4]^2 / 2)
grad_Linf_m <- c(0,0, Linf_m, par_m[4] * Linf_m)
var_Linf_m <- t(grad_Linf_m) %*% cov_m %*% grad_Linf_m
se_Linf_m <- sqrt(var_Linf_m)
CI_Linf_m <- Linf_m + c(-1.96, 1.96) * se_Linf_m


############################################################
#  SUMMARY TABLE
############################################################

table_bl <- data.frame(
  Sex        = c("Female", "Male"),
  k          = c(par_f[1], par_m[1]),
  zeta       = c(par_f[2], par_m[2]),
  mu         = c(par_f[3], par_m[3]),
  sigma      = c(par_f[4], par_m[4]),
  AIC        = c(AIC_f, AIC_m),
  E_Linf     = c(Linf_f, Linf_m),
  E_Linf_L95 = c(CI_Linf_f[1], CI_Linf_m[1]),
  E_Linf_U95 = c(CI_Linf_f[2], CI_Linf_m[2]),
  k_L95      = c(CI_k_f[1], CI_k_m[1]),
  k_U95      = c(CI_k_f[2], CI_k_m[2])
)

table_bl[num_cols] <- lapply(table_bl[num_cols], function(x) round(x, 4))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(table_bl, "results/tables/table2_bl_model.csv", row.names = FALSE)
