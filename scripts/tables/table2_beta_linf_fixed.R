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

res_f <- optim(c(200, 0.1, 50), LL1)$par
Linf_f <- res[1]
k_f    <- res[2]
zeta_f <- res[3]

alpha <- (1 - exp(-k*T)) * (zeta-1)
beta  <- exp(-k*T) * (zeta-1)
res_f

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

res_m <- optim(c(200, 0.1, 50), LL1)$par
Linf_m <- res[1]
k_m    <- res[2]
zeta_m <- res[3]

alpha <- (1 - exp(-k*T)) * (zeta-1)
beta  <- exp(-k*T) * (zeta-1)
res_m

# --- Female estimates ---
res_f <- optim(c(200, 0.1, 50), LL1)$par
Linf_f <- res_f[1]; k_f <- res_f[2]; zeta_f <- res_f[3]

# --- Male estimates ---
res_m <- optim(c(200, 0.1, 50), LL1)$par
Linf_m <- res_m[1]; k_m <- res_m[2]; zeta_m <- res_m[3]

# --- Combine into a results data frame ---
table2_results <- data.frame(
  Sex   = c("Female", "Male"),
  Linf  = c(Linf_f, Linf_m),
  k     = c(k_f, k_m),
  zeta  = c(zeta_f, zeta_m)
)

# --- Save to results/tables ---
write.csv(
  table2_results,
  "results/tables/table2_bg_model.csv",
  row.names = FALSE
)
