dat<-lobster
# ff_GG as you provided
ff_GG <- function(Linf, k, alpha, a_L, b_L, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    mu <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
    shape <- alpha
    scale <- mu / alpha
    log_inc <- dgamma(Xi[j], shape=shape, scale=scale, log=TRUE)
    log_prior <- dgamma(Linf, shape=a_L, rate=b_L, log=TRUE)
    gLinf <- gLinf + log_inc + log_prior
  }
  return(gLinf)
}

# Normalization for integration
gg_GG <- function(Linf, k, alpha, a_L, b_L, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff_GG(Linf, k, alpha, a_L, b_L, Xi, Ti, Li)
  Mi <- optimize(ff_GG, c(MinLinf, MaxLinf), k=k, alpha=alpha, a_L=a_L, b_L=b_L, Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
  return(exp(res - Mi))
}

LL_GG <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  alpha <- theta[2]
  a_L <- theta[3]
  b_L <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    Mi <- optimize(ff_GG, c(MinLinf, MaxLinf), k=k, alpha=alpha, a_L=a_L, b_L=b_L, Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
    res.int <- integrate(gg_GG, lower=MinLinf, upper=MaxLinf, k=k, alpha=alpha, a_L=a_L, b_L=b_L, Xi=Xi, Ti=Ti, Li=Li, MinLinf=MinLinf, MaxLinf=MaxLinf)$value
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

MinLinf <- 160
MaxLinf <- 200
init_theta <- c(0.1, 4, 5, 0.1) # k, alpha, a_L, b_L (set sensible starting values)

dat_female <- lobster[lobster$SEX == 1, ]

res_female_GG <- optim(init_theta, LL_GG, dat=dat_female,
                       MinLinf=MinLinf, MaxLinf=MaxLinf,
                       lower=c(0.01, 0.1, 0.1, 0.01),
                       upper=c(1, 100, 100, 10),
                       method="L-BFGS-B",
                       control=list(maxit=1000, trace=1))

val_female_GG <- res_female_GG$par
neg_logL_female_GG <- res_female_GG$value
logL_female_GG <- -neg_logL_female_GG
num_params <- 4
AIC_female_GG <- 2 * num_params + 2 * neg_logL_female_GG
cat("AIC (GG, Female):", AIC_female_GG, "\n")

dat_male <- lobster[lobster$SEX == 2, ]
res_male_GG <- optim(init_theta, LL_GG, dat=dat_male,
                     MinLinf=MinLinf, MaxLinf=MaxLinf,
                     lower=c(0.01, 0.1, 0.1, 0.01),
                     upper=c(1, 100, 100, 10),
                     method="L-BFGS-B",
                     control=list(maxit=1000, trace=1))
val_male_GG <- res_male_GG$par
neg_logL_male_GG <- res_male_GG$value
logL_male_GG <- -neg_logL_male_GG
AIC_male_GG <- 2 * num_params + 2 * neg_logL_male_GG
cat("AIC (GG, Male):", AIC_male_GG, "\n")
