test_table6ind <- function() {
  data <- lobster
  
  independent_likelihood_kMI <- function(params, data) {
    Linf   <- params[1]; k <- params[2]; lambda <- params[3]
    beta0  <- params[4]; beta1 <- params[5]; beta2 <- params[6]; phi <- params[7]
    loglik_MI <- 0; loglik_IP <- 0
    
    for (i in 1:nrow(data)) {
      Lp <- data$PL[i]; I <- data$INC[i]; T <- data$INT[i]/365.25
      if (T <= 0) next
      if ((Linf - Lp) <= 0) return(1e10)
      mu_I <- (Linf - Lp) * (1 - exp(-k * T))
      if (mu_I <= 0 || !is.finite(mu_I)) return(1e10)
      shape_I <- lambda * mu_I; scale_I <- 1/lambda
      if (shape_I <= 0 || scale_I <= 0) return(1e10)
      loglik_MI <- loglik_MI + dgamma(I, shape=shape_I, scale=scale_I, log=TRUE)
      mu_T <- exp(beta0 + beta1*Lp + beta2*(Lp^2))
      if (!is.finite(mu_T) || mu_T <= 0) return(1e10)
      shape_T <- phi; scale_T <- mu_T/phi
      if (shape_T <= 0 || scale_T <= 0) return(1e10)
      loglik_IP <- loglik_IP + dgamma(T, shape=shape_T, scale=scale_T, log=TRUE)
    }
    negLL <- -(loglik_MI + loglik_IP)
    if (!is.finite(negLL)) return(1e10)
    negLL
  }
  
  # Prepare female/male data
  data_female <- subset(data, SEX == 1)
  data_male   <- subset(data, SEX == 2)
  
  # Initial values and bounds
  init_params <- c(Linf=185, k=0.3, lambda=1, beta0=-1, beta1=0, beta2=0, phi=2)
  lower_bounds <- c(150, 0.01, 0.01, -10, -0.01, -1e-4, 0.1)
  upper_bounds <- c(220, 1.00, 20, 10, 0.01, 1e-4, 20)
  
  # Female fit
  fit_female <- optim(par=init_params, fn=independent_likelihood_kMI, data=data_female,
                      method="L-BFGS-B", lower=lower_bounds, upper=upper_bounds,
                      control=list(maxit=1000))
  AIC_female <- 2*fit_female$value + 2*length(fit_female$par)
  
  # Male fit
  fit_male <- optim(par=init_params, fn=independent_likelihood_kMI, data=data_male,
                    method="L-BFGS-B", lower=lower_bounds, upper=upper_bounds,
                    control=list(maxit=1000))
  AIC_male <- 2*fit_male$value + 2*length(fit_male$par)
  
  # Results table
  results_table <- data.frame(
    Sex      = c("Female","Male"),
    Linf     = c(fit_female$par["Linf"],  fit_male$par["Linf"]),
    k        = c(fit_female$par["k"],     fit_male$par["k"]),
    lambda   = c(fit_female$par["lambda"],fit_male$par["lambda"]),
    beta0    = c(fit_female$par["beta0"], fit_male$par["beta0"]),
    beta1    = c(fit_female$par["beta1"], fit_male$par["beta1"]),
    beta2    = c(fit_female$par["beta2"], fit_male$par["beta2"]),
    phi      = c(fit_female$par["phi"],   fit_male$par["phi"]),
    NegLL    = c(fit_female$value,        fit_male$value),
    AIC      = c(AIC_female,              AIC_male),
    Converged= c(fit_female$convergence==0, fit_male$convergence==0)
  )
  
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive=TRUE)  
  write.csv(results_table, "results/tables/table6_independent_model.csv", row.names=FALSE)
  invisible(results_table)
}
