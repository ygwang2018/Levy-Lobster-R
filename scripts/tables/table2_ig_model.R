f <- dat[dat$SEX == 1, ]      # Sex = 1 (Female); change if needed
T <- f$INT / 365.25           # years
I <- f$INC                    # observed MI
L <- f$PL                     # premoult length

log_dinvgauss <- function(x, mu, lambda) {
  # density:
  # f(x) = sqrt(lambda/(2*pi*x^3)) * exp( -lambda (x-mu)^2 / (2 mu^2 x) )
  0.5 * (log(lambda) - log(2*pi) - 3*log(x)) -
    lambda * (x - mu)^2 / (2 * mu^2 * x)
}

LL_IG <- function(d) {
  
  Linf   <- d[1]
  lambda <- d[2]
  k      <- d[3]
  
  # ---- constraints ----
  if (Linf <= max(L) + 1e-6) return(1e10)
  if (lambda <= 0) return(1e10)
  if (k <= 0) return(1e10)
  
  # mean increment mu_i
  mu <- (Linf - L) * (1 - exp(-k * T))
  if (any(mu <= 0)) return(1e10)
  if (any(I <= 0)) return(1e10)
  
  # log-likelihood
  ll <- log_dinvgauss(I, mu, lambda)
  
  return(-sum(ll))   # minimisation target
}

# Fit IG MI model
res_f <- optim(
  par    = c(180, 1, 0.1),
  fn     = LL_IG,
  method = "L-BFGS-B",
  lower  = c(100, 0.01, 0.01),
  upper  = c(300, 10, 0.4),
  control = list(maxit = 2000)
)

val_f <- res_f$par
negLL_f <- res_f$value
p_f <- length(val_f)
AIC_f <- 2*negLL_f + 2*p_f


########. MALE

dat<-lobster
f <- dat[dat$SEX ==2, ]      # Sex = 2 (Male); change if needed
T <- f$INT / 365.25           # years
I <- f$INC                    # observed MI
L <- f$PL                     # premoult length

log_dinvgauss <- function(x, mu, lambda) {
  # density:
  # f(x) = sqrt(lambda/(2*pi*x^3)) * exp( -lambda (x-mu)^2 / (2 mu^2 x) )
  0.5 * (log(lambda) - log(2*pi) - 3*log(x)) -
    lambda * (x - mu)^2 / (2 * mu^2 * x)
}

LL_IG <- function(d) {
  
  Linf   <- d[1]
  lambda <- d[2]
  k      <- d[3]
  
  # ---- constraints ----
  if (Linf <= max(L) + 1e-6) return(1e10)
  if (lambda <= 0) return(1e10)
  if (k <= 0) return(1e10)
  
  # mean increment mu_i
  mu <- (Linf - L) * (1 - exp(-k * T))
  if (any(mu <= 0)) return(1e10)
  if (any(I <= 0)) return(1e10)
  
  # log-likelihood
  ll <- log_dinvgauss(I, mu, lambda)
  
  return(-sum(ll))   # minimisation target
}

# Fit IG MI model
res_m <- optim(
  par    = c(160, 1, 0.1),
  fn     = LL_IG,
  method = "L-BFGS-B",
  lower  = c(160, 0.01, 0.01),
  upper  = c(180, 10, 1),
  control = list(maxit = 2000)
)

val_m <- res_m$par
negLL_m <- res_m$value
p_m <- length(val_m)
AIC_m <- 2*negLL_m + 2*p_m

# Extract parameter estimates
Linf_IG   <- res_m$par[1]
lambda_IG <- res_m$par[2]
k_IG      <- res_m$par[3]
res_m


# --- Combine into results table ---
table2_results <- data.frame(
  Sex    = c("Female", "Male"),
  Linf   = c(val_f[1], val_m[1]),
  lambda = c(val_f[2], val_m[2]),
  k      = c(val_f[3], val_m[3]),
  AIC    = c(AIC_f, AIC_m)
)

# --- Save to results/tables ---
write.csv(
  table2_results,
  "results/tables/table2_ig_model.csv",
  row.names = FALSE
)
