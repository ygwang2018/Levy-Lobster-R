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
#lambda0 <- 0.1
#k0      <- 0.2
#alpha0  <- 5
#beta0   <- 0.05

#MinLinf <- 190
#MaxLinf <- 230

#val <- optim(
#  c(lambda0, k0, alpha0, beta0),
#  LL,
#  dat=dat,
#  MinLinf=MinLinf,
#  MaxLinf=MaxLinf,
#  method="L-BFGS-B",
#  lower=c(0.01,0.01,0.10,0.01),
#  upper=c(1,0.5,10,100)
#)
#val

dat_f <- subset(lobster, SEX == 1)
MinLinf_f <- max(dat_f$PL) + 0.1
MaxLinf_f <- MinLinf_f + 80

start_f <- c(lambda=0.1, k=0.2, alpha=5, beta=0.05)
fit_f <- optim(start_f, LL, dat=dat_f, MinLinf=MinLinf_f, MaxLinf=MaxLinf_f,
               method="L-BFGS-B", lower=c(0.01,0.01,0.10,0.01), upper=c(1,0.5,10,100),
               control=list(maxit=2000, trace=1))
par_f  <- fit_f$par
negLL_f <- fit_f$value
AIC_f <- 2*negLL_f + 2*length(par_f)
E_Linf_f <- par_f["alpha"] / par_f["beta"]  # Gamma mean

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

dat_m <- subset(lobster, SEX == 2)
MinLinf_m <- max(dat_m$PL) + 0.1
MaxLinf_m <- MinLinf_m + 80

start_m <- c(lambda=0.1, k=0.2, alpha=5, beta=0.05)
fit_m <- optim(start_m, LL, dat=dat_m, MinLinf=MinLinf_m, MaxLinf=MaxLinf_m,
               method="L-BFGS-B", lower=c(0.01,0.01,0.10,0.01), upper=c(1,0.5,10,100),
               control=list(maxit=2000, trace=1))
par_m  <- fit_m$par
negLL_m <- fit_m$value
AIC_m <- 2*negLL_m + 2*length(par_m)
E_Linf_m <- par_m["alpha"] / par_m["beta"]


#val <- optim(
 # c(lambda0, k0, alpha0, beta0),
# LL,
# dat=dat,
  #MinLinf=MinLinf,
  #MaxLinf=MaxLinf,
  #method="L-BFGS-B",
  #lower=c(0.01,0.01,0.10,0.01),
  #upper=c(1,0.5,10,100)
#)
#val

# Helper functions (your originals)
log_MI_gamma <- function(Delta, Linf, lambda, k, Li, Ti) {
  Iij <- (Linf - Li) * (1 - exp(-k * Ti))
  shape <- lambda * Iij
  rate  <- lambda
  if (Iij <= 0 || Delta <= 0 || shape <= 0 || rate <= 0) return(-Inf)
  shape*log(rate) + (shape-1)*log(Delta) - rate*Delta - lgamma(shape)
}
log_prior_Linf <- function(Linf, alpha, beta) dgamma(Linf, shape=alpha, rate=beta, log=TRUE)
loglik_per_individual <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta) {
  s <- log_prior_Linf(Linf, alpha, beta)
  for (j in seq_along(Xi)) s <- s + log_MI_gamma(Xi[j], Linf, lambda, k, Li[j], Ti[j])
  s
}
integrand <- function(Linf, Xi, Li, Ti, lambda, k, alpha, beta) {
  val <- loglik_per_individual(Linf, Xi, Li, Ti, lambda, k, alpha, beta)
  exp(val)
}
LL <- function(theta, dat, MinLinf, MaxLinf) {
  lambda <- theta[1]; k <- theta[2]; alpha <- theta[3]; beta <- theta[4]
  if (lambda <= 0 || k <= 0 || alpha <= 0 || beta <= 0) return(1e12)
  LL_total <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT / 365.25
    # Optional skip: increments exceeding any Linf bound
    if (any(Li + Xi >= MaxLinf)) next
    val <- try(integrate(integrand, lower=MinLinf, upper=MaxLinf,
                         Xi=Xi, Li=Li, Ti=Ti, lambda=lambda, k=k, alpha=alpha, beta=beta),
               silent=TRUE)
    if (inherits(val, "try-error") || val$value <= 0 || !is.finite(val$value)) return(1e12)
    LL_total <- LL_total + log(val$value)
  }
  -LL_total
}


# --- Combine and save ---
table_results <- data.frame(
  Sex    = c("Female","Male"),
  lambda = c(par_f["lambda"], par_m["lambda"]),
  k      = c(par_f["k"],      par_m["k"]),
  alpha  = c(par_f["alpha"],  par_m["alpha"]),
  beta   = c(par_f["beta"],   par_m["beta"]),
  E_Linf = c(E_Linf_f,        E_Linf_m),
  AIC    = c(AIC_f,           AIC_m)
)

# Optional rounding for presentation
table_results_round <- transform(
  table_results,
  lambda = round(lambda, 4),
  k      = round(k, 4),
  alpha  = round(alpha, 3),
  beta   = round(beta, 3),
  E_Linf = round(E_Linf, 2),
  AIC    = round(AIC, 1)
)

write.csv(table_results_round, "results/tables/table2_igg_model.csv", row.names = FALSE)
