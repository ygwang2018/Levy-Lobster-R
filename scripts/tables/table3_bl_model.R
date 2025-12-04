#    FEMALE BL MODEL
dat<-lobster

# SEX == 1 (Female)
dat <- subset(dat, SEX == 1)
MinLinf <- max(dat$PL + dat$INC) + 0.1   
MaxLinf <- MinLinf + 80                  
ff_ind <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
  
  if (!is.finite(Linf) || Linf <= max(Li) + 1e-8)
    return(-1e10)   # soft penalty
  
  logPrior <- -log(Linf) - log(sigma) - 0.5*log(2*pi) -
    ( (log(Linf) - mu)^2 / (2*sigma^2) )
  
  out <- logPrior
  
  for (j in seq_along(Xi)) {
    
    denom <- Linf - Li[j]
    if (denom <= 1e-8)
      return(-1e10)
    
    x <- Xi[j] / denom
    if (x <= 1e-8 || x >= 1 - 1e-8)
      return(-1e10)
    
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    
    if (a <= 1e-8 || b <= 1e-8)
      return(-1e10)
    
    out <- out +
      dbeta(x, a, b, log = TRUE) -
      log(denom)
  }
  
  return(out)
}
kernel_norm <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
  val <- ff_ind(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
  exp(val - Mi)   # ALWAYS finite
}
max_over_Linf <- function(fun, lower, upper) {
  vals <- c(
    optimize(fun, interval = c(lower, upper), maximum = TRUE)$objective,
    fun(lower),
    fun(upper),
    fun((lower + upper)/2)
  )
  max(vals)
}
safe_integrate <- function(fun, lower, upper, n = 100) {
  xs <- seq(lower, upper, length.out = n)
  vals <- sapply(xs, fun)
  dx <- (upper - lower) / (n - 1)
  sum(vals) * dx
}
LL <- function(theta, dat, MinLinf, MaxLinf) {
  
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
    
    Mi <- max_over_Linf(
      function(Ls) ff_ind(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    val <- safe_integrate(
      function(Ls)
        kernel_norm(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    if (val <= 0 || !is.finite(val))
      return(1e12)
    
    LLtot <- LLtot + log(val) + Mi
  }
  
  return(-LLtot)
}
init <- c(k = 0.10, zeta = 6, mu = log(180), sigma = 0.20)

fit <- optim(
  par     = init,
  fn      = LL,
  dat     = dat,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.1, 1.01, log(120), 0.05),
  upper   = c(0.5,200, log(260), 1.0),
  control = list(maxit = 300)
)

fit$par
fit$value
negLL <- fit$value
p     <- length(init)

AIC_value <- 2*negLL + 2*p
AIC_value
sig<- fit$par[4]
Linf<-exp(fit$par[3]+sig^2/2)
Linf

####FEMALE CI#######

if (!exists("fit")) stop("fit (female model) not found")

# Extract MLEs
par_f <- fit$par
k_f     <- par_f[1]
zeta_f  <- par_f[2]
mu_f    <- par_f[3]
sigma_f <- par_f[4]

H_f <- optimHess(
  par = par_f,
  fn  = LL,
  dat = dat,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf
)

cov_f <- solve(H_f)
se_f  <- sqrt(diag(cov_f))

# 95% CI for k
CI_k_f <- c(
  k_f - 1.96 * se_f[1],
  k_f + 1.96 * se_f[1]
)
names(CI_k_f) <- c("Lower 95%", "Upper 95%")
CI_k_f

# 95% CI for E[Linf] 

Linf_f <- exp(mu_f + sigma_f^2 / 2)

grad_Linf_f <- c(0, 0,
                 Linf_f,
                 sigma_f * Linf_f)

var_Linf_f <- t(grad_Linf_f) %*% cov_f %*% grad_Linf_f
se_Linf_f  <- sqrt(var_Linf_f)

CI_Linf_f <- c(
  Linf_f - 1.96 * se_Linf_f,
  Linf_f + 1.96 * se_Linf_f
)
names(CI_Linf_f) <- c("Lower 95%", "Upper 95%")

CI_Linf_f

###########################
#    MALE BL MODEL

dat<-lobster
dat <- dat[dat$SEX == 2, ]
MinLinf <- max(dat$PL + dat$INC) + 0.1   
MaxLinf <- MinLinf + 80                  
ff_ind <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
  
  if (!is.finite(Linf) || Linf <= max(Li) + 1e-8)
    return(-1e10)   # soft penalty
  
  logPrior <- -log(Linf) - log(sigma) - 0.5*log(2*pi) -
    ( (log(Linf) - mu)^2 / (2*sigma^2) )
  
  out <- logPrior
  
  for (j in seq_along(Xi)) {
    
    denom <- Linf - Li[j]
    if (denom <= 1e-8)
      return(-1e10)
    
    x <- Xi[j] / denom
    if (x <= 1e-8 || x >= 1 - 1e-8)
      return(-1e10)
    
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    
    if (a <= 1e-8 || b <= 1e-8)
      return(-1e10)
    
    out <- out +
      dbeta(x, a, b, log = TRUE) -
      log(denom)
  }
  
  return(out)
}
kernel_norm <- function(Linf, Mi, k, zeta, mu, sigma, Xi, Ti, Li) {
  val <- ff_ind(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
  exp(val - Mi)   # ALWAYS finite
}
max_over_Linf <- function(fun, lower, upper) {
  vals <- c(
    optimize(fun, interval = c(lower, upper), maximum = TRUE)$objective,
    fun(lower),
    fun(upper),
    fun((lower + upper)/2)
  )
  max(vals)
}
safe_integrate <- function(fun, lower, upper, n = 100) {
  xs <- seq(lower, upper, length.out = n)
  vals <- sapply(xs, fun)
  dx <- (upper - lower) / (n - 1)
  sum(vals) * dx
}
LL <- function(theta, dat, MinLinf, MaxLinf) {
  
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
    
    Mi <- max_over_Linf(
      function(Ls) ff_ind(Ls, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    val <- safe_integrate(
      function(Ls)
        kernel_norm(Ls, Mi, k, zeta, mu, sigma, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    if (val <= 0 || !is.finite(val))
      return(1e12)
    
    LLtot <- LLtot + log(val) + Mi
  }
  
  return(-LLtot)
}
init <- c(k = 0.010, zeta = 6, mu = log(180), sigma = 0.20)

fit <- optim(
  par     = init,
  fn      = LL,
  dat     = dat,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.01, 1.01, log(120), 0.05),
  upper   = c(0.5,   200, log(260), 1.0),
  control = list(maxit = 300)
)

fit$par
fit$value
negLL <- fit$value
p     <- length(init)

AIC_value <- 2*negLL + 2*p
AIC_value
sig<- fit$par[4]
Linf<-exp(fit$par[3]+sig^2/2)
Linf

####### MALE CI##############

if (!exists("fit")) stop("fit (male model) not found")

par_m <- fit$par
k_m     <- par_m[1]
zeta_m  <- par_m[2]
mu_m    <- par_m[3]
sigma_m <- par_m[4]

# Hessian
H_m <- optimHess(
  par = par_m,
  fn  = LL,
  dat = dat,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf
)

cov_m <- solve(H_m)
se_m  <- sqrt(diag(cov_m))

# 95% CI for k
CI_k_m <- c(
  k_m - 1.96 * se_m[1],
  k_m + 1.96 * se_m[1]
)
names(CI_k_m) <- c("Lower 95%", "Upper 95%")
CI_k_m

# CI for E[Linf] 

Linf_m <- exp(mu_m + sigma_m^2 / 2)

grad_Linf_m <- c(0, 0,
                 Linf_m,
                 sigma_m * Linf_m)

var_Linf_m <- t(grad_Linf_m) %*% cov_m %*% grad_Linf_m
se_Linf_m  <- sqrt(var_Linf_m)

CI_Linf_m <- c(
  Linf_m - 1.96 * se_Linf_m,
  Linf_m + 1.96 * se_Linf_m
)
names(CI_Linf_m) <- c("Lower 95%", "Upper 95%")

CI_Linf_m
