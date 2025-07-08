dat<-lobster
dat <- dat[dat$SEX == 1, ]

ff <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li) {
  Linf <- as.numeric(Linf[1])  # force scalar
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    denom <- Linf - Li[j]
    if (!is.finite(denom) || denom <= 0) return(-Inf)
    
    x <- Xi[j] / denom
    if (!is.finite(x) || x <= 0 || x >= 1) return(-Inf)
    
    bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x)
    bb <- bb - lgamma(a) - lgamma(b) - log(sigma) - 0.5 * log(2 * pi)
    bb <- bb - ((log(Linf) - miu)^2) / (2 * sigma^2) - log(denom)
    
    if (!is.finite(bb) || length(bb) != 1) return(-Inf)
    
    gLinf <- gLinf + bb
  }
  return(gLinf)
}

gg <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  Linf <- as.numeric(Linf[1])  # force scalar
  res <- ff(Linf, k, zeta, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(
    function(Linf_scalar) ff(Linf_scalar, k, zeta, miu, sigma, Xi, Ti, Li),
    c(MinLinf, MaxLinf),
    maximum = TRUE
  )$objective
  
  val <- exp(res - Mi)
  if (!is.finite(val) || length(val) != 1) val <- 0
  
  return(as.numeric(val))
}

LL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  zeta <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    Mi <- optimize(
      function(Linf_scalar) ff(Linf_scalar, k, zeta, miu, sigma, Xi, Ti, Li),
      c(MinLinf, MaxLinf),
      maximum = TRUE
    )$objective
    
    res.int <- tryCatch({
      suppressWarnings(
        integrate(
          function(Linf) gg(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf),
          lower = MinLinf,
          upper = MaxLinf,
          subdivisions = 200L,
          abs.tol = 1e-6
        )$value
      )
    }, error = function(e) {
      warning(paste("Integration failed:", e$message))
      return(1e-10)
    })
    
    all.LL <- all.LL + log(res.int) + Mi
  }
  print(c(k, zeta, miu, sigma))
  return(-all.LL)
}

# Initial parameter guesses
k <- 0.1
zeta <- 4
miu <- 1
sigma <- 2.88
MinLinf <- 160
MaxLinf <- 200

# Optimization
val <- optim(c(k, zeta, miu, sigma), LL, dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf,
             lower = c(0.01, 0.1, 0.1, 0.1), upper = c(1, 100, 10, 10),
             method = "L-BFGS-B")$par

# Final parameter values
val

# Derived Linf estimate
k <- val[1]
zeta <- val[2]
miu <- val[3]
sigma <- val[4]
Linf <- exp(miu + (sigma^2) / 2)
Linf


#————————————————————————————————————————————————————————

#male
dat<-lobster
dat <- dat[dat$SEX == 2, ]
# Filter to male lobsters
dat <- dat[dat$SEX == 2, ]

ff <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li) {
  Linf <- as.numeric(Linf[1])  # force scalar
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    denom <- Linf - Li[j]
    if (!is.finite(denom) || denom <= 0) return(-Inf)
    
    x <- Xi[j] / denom
    if (!is.finite(x) || x <= 0 || x >= 1) return(-Inf)
    
    lnx <- log(x)
    ln1mx <- log(1 - x)
    
    bb <- lgamma(a + b) + (a - 1) * lnx + (b - 1) * ln1mx
    bb <- bb - lgamma(a) - lgamma(b) - log(sigma) - 0.5 * log(2 * pi)
    bb <- bb - ((log(Linf) - miu)^2) / (2 * sigma^2) - log(denom)
    
    if (!is.finite(bb) || length(bb) != 1) return(-Inf)
    
    gLinf <- gLinf + bb
  }
  return(gLinf)
}

gg <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  Linf <- as.numeric(Linf[1])  # force scalar
  res <- ff(Linf, k, zeta, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(
    function(Linf_scalar) ff(Linf_scalar, k, zeta, miu, sigma, Xi, Ti, Li),
    c(MinLinf, MaxLinf),
    maximum = TRUE
  )$objective
  
  val <- exp(res - Mi)
  if (!is.finite(val) || length(val) != 1) val <- 0
  
  return(as.numeric(val))
}

LL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  zeta <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    Mi <- optimize(
      function(Linf_scalar) ff(Linf_scalar, k, zeta, miu, sigma, Xi, Ti, Li),
      c(MinLinf, MaxLinf),
      maximum = TRUE
    )$objective
    
    res.int <- tryCatch({
      suppressWarnings(
        integrate(
          function(Linf) gg(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf),
          lower = MinLinf,
          upper = MaxLinf,
          subdivisions = 200L,
          abs.tol = 1e-6
        )$value
      )
    }, error = function(e) {
      warning(paste("Integration failed:", e$message))
      return(1e-10)
    })
    
    all.LL <- all.LL + log(res.int) + Mi
  }
  print(c(k, zeta, miu, sigma))
  return(-all.LL)
}

# Initial parameter guesses
k <- 0.1
zeta <- 4
miu <- 1
sigma <- 2.88
MinLinf <- max(dat$PL) + 1
MaxLinf <- MinLinf + 40

# Optimization
val <- optim(c(k, zeta, miu, sigma), LL, dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf,
             lower = c(0.01, 0.1, 0.1, 0.1), upper = c(1, 100, 10, 10),
             method = "L-BFGS-B")$par

# Final parameter values
val

# Derived Linf estimate
k <- val[1]
zeta <- val[2]
miu <- val[3]
sigma <- val[4]
Linf <- exp(miu + (sigma^2) / 2)
Linf
