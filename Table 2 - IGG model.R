dat<-lobster
dat <- dat[dat$SEX == 1, ]
#dat <- dat[order(dat$CL),]

#  Gamma prior for Linf ~ Gamma(alpha, beta)
#  Gamma MI model:  Δ | Linf ~ Gamma(shape = λ * Iij, rate=λ)
#  Iij = (Linf - Li) * (1 - exp(-k*T))


# MI log-density for ONE observation j
log_MI_gamma <- function(Delta, Linf, lambda, k, Li, Ti)
{
  Iij <- (Linf - Li) * (1 - exp(-k * Ti))
  shape <- lambda * Iij
  rate  <- lambda
  
  # log f(Δ_ij | Linf)
  out <- shape*log(rate) + (shape-1)*log(Delta) - rate*Delta - lgamma(shape)
  return(out)
}

# log prior for ONE Linf
log_prior_Linf <- function(Linf, alpha, beta)
{
  dgamma(Linf, shape=alpha, rate=beta, log=TRUE)
}

# full log-likelihood for ONE individual i
loglik_per_individual <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
{
  s <- 0
  n <- length(Xi)
  
  for (j in 1:n)
    s <- s + log_MI_gamma(Xi[j], Linf, lambda, k, Li[j], Ti[j])
  
  # add PRIOR only once per individual
  s <- s + log_prior_Linf(Linf, alpha, beta)
  return(s)
}

# to integrate over Linf ∈ [MinLinf, MaxLinf]:
integrand <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
{
  exp(loglik_per_individual(Linf, Xi, Li, Ti, lambda, k, alpha, beta))
}

#  Full negative log-likelihood

LL <- function(theta, dat, MinLinf, MaxLinf)
{
  lambda <- theta[1]
  k      <- theta[2]
  alpha  <- theta[3]
  beta   <- theta[4]
  
  individuals <- unique(dat$LOBSTER)
  LL_total <- 0
  
  for (iid in individuals)
  {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    # integrate f(Δ, Linf)
    val <- try(
      integrate(integrand,
                lower=MinLinf, upper=MaxLinf,
                Xi=Xi, Li=Li, Ti=Ti,
                lambda=lambda, k=k, alpha=alpha, beta=beta),
      silent=TRUE)
    
    if (class(val) == "try-error" || val$value <= 0)
      return(1e12)
    
    LL_total <- LL_total + log(val$value)
  }
  
  return(-LL_total)
}

# OPTIMIZATION
lambda0 <- 0.1
k0      <- 0.2
alpha0  <- 5
beta0   <- 0.05

MinLinf <- 190
MaxLinf <- 230

val <- optim(
  c(lambda0, k0, alpha0, beta0),
  LL,
  dat=dat,
  MinLinf=MinLinf,
  MaxLinf=MaxLinf,
  method="L-BFGS-B",
  lower=c(0.01,0.01,0.10,0.01),
  upper=c(1,0.5,10,100)
)

val

####### MALE

dat <- dat[dat$SEX == 2, ]
dat <- dat[order(dat$PL),]

# MI log-density for ONE observation j
log_MI_gamma <- function(Delta, Linf, lambda, k, Li, Ti)
{
  Iij <- (Linf - Li) * (1 - exp(-k * Ti))
  shape <- lambda * Iij
  rate  <- lambda
  
  # log f(Δ_ij | Linf)
  out <- shape*log(rate) + (shape-1)*log(Delta) - rate*Delta - lgamma(shape)
  return(out)
}

# log prior for ONE Linf
log_prior_Linf <- function(Linf, alpha, beta)
{
  dgamma(Linf, shape=alpha, rate=beta, log=TRUE)
}

# full log-likelihood for ONE individual i
loglik_per_individual <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
{
  s <- 0
  n <- length(Xi)
  
  for (j in 1:n)
    s <- s + log_MI_gamma(Xi[j], Linf, lambda, k, Li[j], Ti[j])
  
  # add PRIOR only once per individual
  s <- s + log_prior_Linf(Linf, alpha, beta)
  return(s)
}

# to integrate over Linf ∈ [MinLinf, MaxLinf]:
integrand <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
{
  exp(loglik_per_individual(Linf, Xi, Li, Ti, lambda, k, alpha, beta))
}

#  Full negative log-likelihood
LL <- function(theta, dat, MinLinf, MaxLinf)
{
  lambda <- theta[1]
  k      <- theta[2]
  alpha  <- theta[3]
  beta   <- theta[4]
  
  individuals <- unique(dat$LOBSTER)
  LL_total <- 0
  
  for (iid in individuals)
  {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    # integrate f(Δ, Linf)
    val <- try(
      integrate(integrand,
                lower=MinLinf, upper=MaxLinf,
                Xi=Xi, Li=Li, Ti=Ti,
                lambda=lambda, k=k, alpha=alpha, beta=beta),
      silent=TRUE)
    
    if (class(val) == "try-error" || val$value <= 0)
      return(1e12)
    
    LL_total <- LL_total + log(val$value)
  }
  
  return(-LL_total)
}

# OPTIMIZATION

lambda0 <- 0.1
k0      <- 0.2
alpha0  <- 5
beta0   <- 0.05

MinLinf <- 190
MaxLinf <- 230

val <- optim(
  c(lambda0, k0, alpha0, beta0),
  LL,
  dat=dat,
  MinLinf=MinLinf,
  MaxLinf=MaxLinf,
  method="L-BFGS-B",
  lower=c(0.01,0.01,0.10,0.01),
  upper=c(1,0.5,10,100)
)

val
