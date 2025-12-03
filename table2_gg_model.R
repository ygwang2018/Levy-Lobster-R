
ff_GG <- function(Linf, k, alpha, a_L, b_L, Xi, Ti, Li) {
  
  # L∞ must be > all premoult lengths
  if (!is.finite(Linf) || Linf <= max(Li) + 1e-10)
    return(-1e12)
  
  # Gamma prior on L∞ (shape=a_L, rate=b_L)
  log_prior <- dgamma(Linf, shape = a_L, rate = b_L, log = TRUE)
  ll <- log_prior
  
  # Loop over increments
  for (j in seq_along(Xi)) {
    
    mu <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))   # mean increment
    if (mu <= 0 || Xi[j] <= 0) return(-1e12)
    
    # Gamma MI: shape = alpha, scale = mu/alpha  (mean = mu)
    shape <- alpha
    scale <- mu / alpha
    
    if (shape <= 0 || scale <= 0) return(-1e12)
    
    log_inc <- dgamma(Xi[j], shape = shape, scale = scale, log = TRUE)
    ll <- ll + log_inc
  }
  
  return(ll)
}

kernel_GG <- function(Linf, Mi, k, alpha, a_L, b_L, Xi, Ti, Li) {
  val <- ff_GG(Linf, k, alpha, a_L, b_L, Xi, Ti, Li)
  exp(val - Mi)
}

safe_integrate <- function(fun, lower, upper, n = 1200) {
  xs <- seq(lower, upper, length.out = n)
  vals <- sapply(xs, fun)
  dx <- (upper - lower) / (n - 1)
  sum(vals) * dx
}

LL_GG <- function(theta, dat, MinLinf, MaxLinf) {
  
  k     <- theta[1]
  alpha <- theta[2]
  a_L   <- theta[3]
  b_L   <- theta[4]
  
  # Basic parameter constraints
  if (k <= 0 || alpha <= 0 || a_L <= 0 || b_L <= 0)
    return(1e12)
  
  LLtot <- 0
  
  for (id in unique(dat$LOBSTER)) {
    
    dati <- subset(dat, LOBSTER == id)
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    # Mode of the integrand for stabilisation
    Mi <- optimize(
      function(Ls) ff_GG(Ls, k, alpha, a_L, b_L, Xi, Ti, Li),
      interval = c(MinLinf, MaxLinf),
      maximum = TRUE
    )$objective
    
    # Integrate normalised kernel
    val <- safe_integrate(
      function(Ls)
        kernel_GG(Ls, Mi, k, alpha, a_L, b_L, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    if (!is.finite(val) || val <= 0)
      return(1e12)
    
    LLtot <- LLtot + log(val) + Mi
  }
  
  return(-LLtot)
}


dat_female <- subset(lobster, SEX == 1)

MinLinf <- max(dat_female$PL) + 0.1
MaxLinf <- MinLinf + 80

init_theta_female <- c(
  k     = 0.15,
  alpha = 10,
  a_L   = 180,    # gamma prior mean = a_L / b_L
  b_L   = 1
)

res_female_GG <- optim(
  par     = init_theta_female,
  fn      = LL_GG,
  dat     = dat_female,
  MinLinf = MinLinf, 
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.01, 0.1, 1e-3, 1e-3),
  upper   = c(1.0, 100,  500, 20),
  control = list(maxit = 2000, trace = 1)
)

val_female_GG <- res_female_GG$par
negLL_female  <- res_female_GG$value
p_female      <- length(val_female_GG)

AIC_female_GG <- 2*negLL_female + 2*p_female
Linf_female_mean <- val_female_GG[3] / val_female_GG[4]   # Gamma mean

#======== FEMALE 
print(val_female_GG)
print("AIC:", AIC_female_GG)
print("E[Linf]:", Linf_female_mean)

# FIT MALE 
dat_male <- subset(lobster, SEX == 2)

MinLinf <- max(dat_male$PL) + 0.1
MaxLinf <- MinLinf + 80

init_theta_male <- c(
  k     = 0.15,
  alpha = 10,
  a_L   = 160,
  b_L   = 1
)

res_male_GG <- optim(
  par     = init_theta_male,
  fn      = LL_GG,
  dat     = dat_male,
  MinLinf = MinLinf, 
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.01, 0.1, 1e-3, 1e-3),
  upper   = c(1.0, 100,  500, 20),
  control = list(maxit = 2000, trace = 1)
)

val_male_GG <- res_male_GG$par
negLL_male  <- res_male_GG$value
p_male      <- length(val_male_GG)

AIC_male_GG <- 2*negLL_male + 2*p_male
Linf_male_mean <- val_male_GG[3] / val_male_GG[4]   # Gamma mean

###== MALE 
print(val_male_GG)
print("AIC:", AIC_male_GG)
print("E[Linf]:", Linf_male_mean)
