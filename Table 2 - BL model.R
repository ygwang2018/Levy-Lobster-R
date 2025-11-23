dat<-lobster

# SEX == 1 (Female)
dat <- subset(dat, SEX == 1)
MinLinf <- 100    # biologically required
MaxLinf <- 200                   # wide enough
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
safe_integrate <- function(fun, lower, upper, n = 200) {
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

#————————————————————————————————————————————————————————

#male

dat<-lobster
dat <- dat[dat$SEX == 2, ]
MinLinf <- max(dat$PL + dat$INC) + 0.1    # biologically required
MaxLinf <- MinLinf + 80                   # wide enough
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
safe_integrate <- function(fun, lower, upper, n = 500) {
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
