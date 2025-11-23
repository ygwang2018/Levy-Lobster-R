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
res <- optim(
  par    = c(180, 1, 0.1),    # starting values: (Linf, lambda, k)
  fn     = LL_IG,
  method = "L-BFGS-B",
  lower  = c(100, 0.01, 0.01),
  upper  = c(300, 10, 0.4),
  control = list(maxit = 2000, trace = 1)
)

# Extract parameter estimates
Linf_IG   <- res$par[1]
lambda_IG <- res$par[2]
k_IG      <- res$par[3]
res

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
res <- optim(
  par    = c(160, 1, 0.1),    # starting values: (Linf, lambda, k)
  fn     = LL_IG,
  method = "L-BFGS-B",
  lower  = c(160, 0.01, 0.01),
  upper = c(180, 10, 1),
  control = list(maxit = 2000, trace = 1)
)

# Extract parameter estimates
Linf_IG   <- res$par[1]
lambda_IG <- res$par[2]
k_IG      <- res$par[3]
res
