ff_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li) {
  
  # lognormal prior for Linf (added ONCE)
  log_prior <- dlnorm(Linf, meanlog = miu, sdlog = sigma, log = TRUE)
  
  gLinf <- log_prior
  ni <- length(Xi)
  
  for (j in 1:ni) {
    # Gamma increment likelihood
    mu    <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))    # mean increment
    shape <- alpha
    scale <- mu / alpha                                # mean = mu, var = mu^2 / alpha
    
    log_inc <- dgamma(Xi[j], shape = shape, scale = scale, log = TRUE)
    
    gLinf <- gLinf + log_inc
  }
  return(gLinf)
}
gg_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff_GL(Linf, k, alpha, miu, sigma, Xi, Ti, Li)
  Mi  <- optimize(ff_GL, c(MinLinf, MaxLinf),
                  k = k, alpha = alpha, miu = miu, sigma = sigma,
                  Xi = Xi, Ti = Ti, Li = Li,
                  maximum = TRUE)$objective
  exp(res - Mi)
}

LL_GL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  alpha <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    # Skip increments that cannot be explained by any Linf
    if (any(Li + Xi >= MaxLinf)) next
    Mi <- optimize(ff_GL, c(MinLinf, MaxLinf), k=k, alpha=alpha, miu=miu, sigma=sigma, Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
    res.int <- tryCatch({
      integrate(gg_GL, lower=MinLinf, upper=MaxLinf, k=k, alpha=alpha, miu=miu, sigma=sigma,
                Xi=Xi, Ti=Ti, Li=Li, MinLinf=MinLinf, MaxLinf=MaxLinf, rel.tol=1e-4)$value
    }, error=function(e) NA)
    if (is.na(res.int) || res.int <= 0) next
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

MinLinf <- 160
MaxLinf <- 200
init_theta <- c(0.1, 4, 5, 0.5) # k, alpha, miu, sigma (set sensible starting values)
num_params <- 4

####### Female
dat_female <- lobster[lobster$SEX == 1, ]
res_female_GL <- optim(init_theta, LL_GL, dat=dat_female,
                       MinLinf=MinLinf, MaxLinf=MaxLinf,
                       lower=c(0.01, 0.1, 0.1, 0.01),
                       upper=c(1, 100, 10, 10),
                       method="L-BFGS-B",
                       control=list(maxit=1000, trace=1))
val_female <- res_female_GL$par
negLL_f    <- res_female_GL$value
AIC_f      <- 2 * num_params + 2 * negLL_f
Linf_f     <- exp(val_female[3] + (val_female[4]^2) / 2)

AIC_f                        
Linf_f                  
val_female
val_female_GL[1]

print("Convergence:", res_female_GL$convergence == 0)
print("AIC:", AIC_female_GL)
print("Estimated Linf:", Linf_female_GL)
print("Parameters: k =", val_female_GL[1], ", alpha =", val_female_GL[2], 
    ", miu =", val_female_GL[3], ", sigma =", val_female_GL[4])

####### Male

dat_male <- lobster[lobster$SEX == 2, ]
res_male_GL <- optim(init_theta, LL_GL, dat=dat_male,
                     MinLinf=MinLinf, MaxLinf=MaxLinf,
                     lower=c(0.01, 0.1, 0.1, 0.01),
                     upper=c(1, 100, 10, 10),
                     method="L-BFGS-B",
                     control=list(maxit=1000, trace=1))
val_male   <- res_male_GL$par
negLL_m    <- res_male_GL$value
AIC_m      <- 2 * num_params + 2 * negLL_m
Linf_m     <- exp(val_male[3] + (val_male[4]^2) / 2)
AIC_m
Linf_m
val_male[1]

print(" Male Results")
print("Convergence:", res_male_GL$convergence == 0)
print("AIC:", AIC_male_GL)
print("Estimated Linf:", Linf_male_GL)
print("Parameters: k =", val_male_GL[1], ", alpha =", val_male_GL[2], 
    ", miu =", val_male_GL[3], ", sigma =", val_male_GL[4])
                        

# --- Combine into results table ---
table2_results <- data.frame(
  Sex    = c("Female", "Male"),
  k      = c(val_female[1], val_male[1]),
  alpha  = c(val_female[2], val_male[2]),
  miu    = c(val_female[3], val_male[3]),
  sigma  = c(val_female[4], val_male[4]),
  E_Linf = c(Linf_f, Linf_m),
  AIC    = c(AIC_f, AIC_m)
)

# --- Save to results/tables ---
write.csv(
  table2_results,
  "results/tables/table2_gl_model.csv",
  row.names = FALSE
)

