####### FEMALE
dat<-lobster
f <- dat[dat$SEX == 1, ]
T <- f$INT / 365.25
I <- f$INC
L <- f$PL
N <- length(f[,1])

LL1 <- function(d)
{
  Linf <- d[1]
  k    <- d[2]
  zeta <- d[3]
  
  # constraints
  if (Linf <= max(L) + 1e-6) return(1e10)
  if (k <= 0) return(1e10)
  if (zeta <= 1) return(1e10)
  
  BG <- 0
  for (i in 1:N) {
    
    # check support
    if (I[i] >= Linf - L[i]) return(1e10)
    
    alpha <- (1 - exp(-k*T[i])) * (zeta-1)
    beta  <- exp(-k*T[i]) * (zeta-1)
    
    # alpha, beta must be >0
    if (alpha <= 0 || beta <= 0) return(1e10)
    
    x <- I[i] / (Linf - L[i])
    
    BG <- BG +
      log(gamma(alpha+beta)) -
      log(gamma(alpha)) - log(gamma(beta)) +
      (alpha-1)*log(x) + (beta-1)*log(1-x) -
      log(Linf - L[i])
  }
  -BG
}

res <- optim(c(200, 0.1, 50), LL1)$par
Linf <- res[1]
k    <- res[2]
zeta <- res[3]

alpha <- (1 - exp(-k*T)) * (zeta-1)
beta  <- exp(-k*T) * (zeta-1)
res

#######.  MALE

####### FEMALE
dat<-lobster
f <- dat[dat$SEX == 2, ]
T <- f$INT / 365.25
I <- f$INC
L <- f$PL
N <- length(f[,1])

LL1 <- function(d)
{
  Linf <- d[1]
  k    <- d[2]
  zeta <- d[3]
  
  # constraints
  if (Linf <= max(L) + 1e-6) return(1e10)
  if (k <= 0) return(1e10)
  if (zeta <= 1) return(1e10)
  
  BG <- 0
  for (i in 1:N) {
    
    # check support
    if (I[i] >= Linf - L[i]) return(1e10)
    
    alpha <- (1 - exp(-k*T[i])) * (zeta-1)
    beta  <- exp(-k*T[i]) * (zeta-1)
    
    # alpha, beta must be >0
    if (alpha <= 0 || beta <= 0) return(1e10)
    
    x <- I[i] / (Linf - L[i])
    
    BG <- BG +
      log(gamma(alpha+beta)) -
      log(gamma(alpha)) - log(gamma(beta)) +
      (alpha-1)*log(x) + (beta-1)*log(1-x) -
      log(Linf - L[i])
  }
  -BG
}

res <- optim(c(200, 0.1, 50), LL1)$par
Linf <- res[1]
k    <- res[2]
zeta <- res[3]

alpha <- (1 - exp(-k*T)) * (zeta-1)
beta  <- exp(-k*T) * (zeta-1)
res
