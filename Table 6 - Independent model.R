# ----------------------------------------
# New independent likelihood function with k in MI
# ----------------------------------------
independent_likelihood_kMI <- function(params, data) {
  Linf  <- params[1]
  k     <- params[2]
  lambda<- params[3]
  beta0 <- params[4]
  beta1 <- params[5]
  beta2 <- params[6]
  phi   <- params[7]
  
  loglik_MI <- 0
  loglik_IP <- 0
  
  for (i in 1:nrow(data)) {
    Lp <- data$PL[i]
    I  <- data$INC[i]
    T  <- data$INT[i] / 365.25
    
    if ((Linf - Lp) <= 0 || T <= 0) return(1e10)
    
    mu_I <- (Linf - Lp) * (1 - exp(-k * T))
    if (mu_I <= 0) return(1e10)
    
    shape_I <- lambda * mu_I
    scale_I <- 1 / lambda
    if (!is.finite(shape_I) || shape_I <= 0 || !is.finite(scale_I) || scale_I <= 0) return(1e10)
    
    log_MI <- dgamma(I, shape = shape_I, scale = scale_I, log = TRUE)
    loglik_MI <- loglik_MI + log_MI
    
    mu_T <- exp(beta0 + beta1 * Lp + beta2 * Lp^2)
    if (!is.finite(mu_T) || mu_T <= 0 || mu_T > 500) return(1e10)
    
    shape_T <- phi
    scale_T <- mu_T / phi
    if (!is.finite(shape_T) || shape_T <= 0 || !is.finite(scale_T) || scale_T <= 0) return(1e10)
    
    log_IP <- dgamma(T, shape = shape_T, scale = scale_T, log = TRUE)
    loglik_IP <- loglik_IP + log_IP
  }
  
  neg_loglik <- -(loglik_MI + loglik_IP)
  if (!is.finite(neg_loglik)) return(1e10)
  
  return(neg_loglik)
}

# ----------------------------------------
# Prepare data subsets
# ----------------------------------------
data_female <- subset(data, SEX == 1)
data_male   <- subset(data, SEX == 2)

# ----------------------------------------
# Starting parameter values (safe)
# ----------------------------------------
init_params <- c(Linf = 185, k = 0.3, lambda = 1, 
                 beta0 = -1, beta1 = 0, beta2 = 0, phi = 2)

lower_bounds <- c(150, 0.01, 0.01, -10, -0.01, -0.0001, 0.1)
upper_bounds <- c(220, 1,   20,    10,  0.01,  0.0001, 20)

# ----------------------------------------
# Female
# ----------------------------------------
cat("\n========== Female ==========\n")

fit_female <- optim(par = init_params,
                    fn = independent_likelihood_kMI,
                    data = data_female,
                    method = "L-BFGS-B",
                    lower = lower_bounds,
                    upper = upper_bounds,
                    control = list(maxit = 1000, trace = 1))

num_params <- length(fit_female$par)
AIC_female <- 2 * fit_female$value + 2 * num_params

cat("AIC    =", AIC_female, "\n")
cat("Neg Log-likelihood =", fit_female$value, "\n")
cat("Convergence =", fit_female$convergence == 0, "\n")
cat("Message:", fit_female$message, "\n")
cat("Linf =", fit_female$par["Linf"], "\n")
cat("k (MI) =", fit_female$par["k"], "\n")
cat("lambda =", fit_female$par["lambda"], "\n")
cat("beta0 =", fit_female$par["beta0"], "\n")
cat("beta1 =", fit_female$par["beta1"], "\n")
cat("beta2 =", fit_female$par["beta2"], "\n")
cat("phi   =", fit_female$par["phi"], "\n")

# ----------------------------------------
# Male
# ----------------------------------------
cat("\n========== Male ==========\n")

fit_male <- optim(par = init_params,
                  fn = independent_likelihood_kMI,
                  data = data_male,
                  method = "L-BFGS-B",
                  lower = lower_bounds,
                  upper = upper_bounds,
                  control = list(maxit = 1000, trace = 1))

AIC_male <- 2 * fit_male$value + 2 * num_params

cat("AIC    =", AIC_male, "\n")
cat("Neg Log-likelihood =", fit_male$value, "\n")
cat("Convergence =", fit_male$convergence == 0, "\n")
cat("Message:", fit_male$message, "\n")
cat("Linf =", fit_male$par["Linf"], "\n")
cat("k (MI) =", fit_male$par["k"], "\n")
cat("lambda =", fit_male$par["lambda"], "\n")
cat("beta0 =", fit_male$par["beta0"], "\n")
cat("beta1 =", fit_male$par["beta1"], "\n")
cat("beta2 =", fit_male$par["beta2"], "\n")
cat("phi   =", fit_male$par["phi"], "\n")
