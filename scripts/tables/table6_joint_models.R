dat<-lobster
ff <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    x <- Xi[j] / (Linf - Li[j])
    bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) - lgamma(a) - lgamma(b) - 
      log(sigma) - 0.5 * log(2 * pi) - ((log(Linf) - miu)^2) / (2 * sigma^2) - log(Linf - Li[j])
    gLinf <- gLinf + bb
  }
  return(gLinf)
}

gg <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff(Linf, k, zeta, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(ff, c(MinLinf, MaxLinf), k = k, zeta = zeta, miu = miu, sigma = sigma,
                 Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
  return(exp(res - Mi))
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
    Mi <- optimize(ff, c(MinLinf, MaxLinf), k = k, zeta = zeta, miu = miu, sigma = sigma,
                   Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
    res.int <- integrate(gg, lower = MinLinf, upper = MaxLinf, k = k, zeta = zeta, miu = miu,
                         sigma = sigma, Xi = Xi, Ti = Ti, Li = Li,
                         MinLinf = MinLinf, MaxLinf = MaxLinf)$value
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

MinLinf <- 160
MaxLinf <- 200
init_theta <- c(0.1, 4, 1, 2.88)

# Female
dat <- lobster
dat_female <- dat[dat$SEX == 1, ]

res_female <- optim(init_theta, LL, dat = dat_female,
                    MinLinf = MinLinf, MaxLinf = MaxLinf,
                    lower = c(0.01, 0.1, 0.1, 0.1),
                    upper = c(1, 100, 10, 10),
                    method = "L-BFGS-B",
                    control = list(maxit = 1000, trace = 1))

val_female <- res_female$par
neg_logL_female <- res_female$value
logL_female <- -neg_logL_female
num_params <- 4
AIC_female <- 2 * (-logL_female) + 2 * num_params
Linf_female <- exp(val_female[3] + (val_female[4]^2) / 2)

#####. Female Results 
print("Convergence:", res_female$convergence == 0)
print( AIC_female)
print( Linf_female)
print(paste( "k =", val_female[1], "zeta =", val_female[2], 
    "miu =", val_female[3], "sigma =", val_female[4]))

# Male
dat <- lobster
dat_male <- dat[dat$SEX == 2, ]

res_male <- optim(init_theta, LL, dat = dat_male,
                  MinLinf = MinLinf, MaxLinf = MaxLinf,
                  lower = c(0.01, 0.1, 0.1, 0.1),
                  upper = c(1, 100, 10, 10),
                  method = "L-BFGS-B",
                  control = list(maxit = 1000, trace = 1))

val_male <- res_male$par
neg_logL_male <- res_male$value
logL_male <- -neg_logL_male 
AIC_male <- 2 * (-logL_male) + 2 * num_params
Linf_male <- exp(val_male[3] + (val_male[4]^2) / 2)

#########. Male Results 
print("Convergence:", res_male$convergence == 0)
print( AIC_male)
print(Linf_male)
print(paste("k =", val_male[1], "zeta =", val_male[2], 
    "miu =", val_male[3],  "sigma =", val_male[4]))

# Assemble results
results_table <- data.frame(
  Sex       = c("Female", "Male"),
  k         = c(val_female[1],   val_male[1]),
  zeta      = c(val_female[2],   val_male[2]),
  miu       = c(val_female[3],   val_male[3]),
  sigma     = c(val_female[4],   val_male[4]),
  E_Linf    = c(Linf_female,     Linf_male),
  LogLik    = c(logL_female,     logL_male),
  NegLogLik = c(neg_logL_female, neg_logL_male),
  AIC       = c(AIC_female,      AIC_male),
  Converged = c(res_female$convergence == 0,
                res_male$convergence == 0)
)

# Save CSV
write.csv(
  results_table,
  file = "results/tables/table6_joint_models.csv",
  row.names = FALSE
)
