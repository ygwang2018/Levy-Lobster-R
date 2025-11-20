ln_params <- function(m, s) {
  sdlog  <- sqrt(log(1 + (s^2 / m^2)))
  meanlog <- log(m) - 0.5 * sdlog^2
  c(meanlog = meanlog, sdlog = sdlog)
}

## Beta–gamma MI–IP simulator for lobster growth
simulate_growth <- function(n, ni_vals, T0_range,
                            Linf_type = c("fixed","lognormal","gamma"),
                            Linf_par, k, zeta, gam, b1, b2) {
  
  Linf_type <- match.arg(Linf_type)
  
  ## number of moults & initial IP
  ni  <- sample(ni_vals, n, replace = TRUE)
  T0  <- sample(T0_range[1]:T0_range[2], n, replace = TRUE) / 365
  
  ## individual L_infty
  Linf_vec <- switch(
    Linf_type,
    "fixed"     = rep(Linf_par, n),
    "lognormal" = rlnorm(n, meanlog = Linf_par[1], sdlog = Linf_par[2]),
    "gamma"     = rgamma(n, shape = Linf_par[1], scale = Linf_par[2])
  )
  
  size_list <- time_list <- vector("list", n)
  
  for (i in 1:n) {
    PL  <- numeric(ni[i] + 1)
    IP  <- numeric(ni[i] + 1)
    INT <- numeric(ni[i] + 1)
    
    PL[1]  <- Linf_vec[i] * (1 - exp(-k * T0[i]))
    IP[1]  <- T0[i]
    INT[1] <- T0[i]
    
    for (j in 1:ni[i]) {
      
      ## beta-subordinator increment using current IP
      shape1 <- (1 - exp(-k * IP[j])) * (zeta - 1)
      shape2 <- exp(-k * IP[j])       * (zeta - 1)
      Lambda <- rbeta(1, shape1, shape2)
      
      PL[j + 1] <- PL[j] + (Linf_vec[i] - PL[j]) * Lambda
      
      ## gamma GLM for next IP
      mu_IP     <- exp(b1 + b2 * PL[j])
      IP[j + 1] <- rgamma(1,
                          shape = gam * mu_IP,
                          scale = 1 / gam)
      INT[j + 1] <- INT[j] + IP[j + 1]
    }
    
    size_list[[i]] <- PL
    time_list[[i]] <- INT
  }
  
  list(time = time_list, size = size_list, Linf = Linf_vec)
}

## Parameters (Table 6 + SD = 12 mm)
set.seed(123)

n      <- 100
years  <- 10
t_grid <- seq(0, years, length.out = 300)

## Females
Linf_f_fixed <- 183.87
f_ln <- ln_params(183.27, 12)
meanlog_f <- f_ln["meanlog"]
sdlog_f   <- f_ln["sdlog"]

k_f    <- 0.283
zeta_f <- 75.055
gam_f  <- 12.535
b1_f   <- -2.087
b2_f   <- 0.011

## Males
Linf_m_fixed <- 228.04
m_ln <- ln_params(184.34, 12)
meanlog_m <- m_ln["meanlog"]
sdlog_m   <- m_ln["sdlog"]

k_m    <- 0.2473
zeta_m <- 59.947
gam_m  <- 11.33
b1_m   <- -1.984
b2_m   <- 0.01

## Panel plotting (NO envelope, NO grey blocks)
plot_panel <- function(res, k, Linf_pop, title) {
  
  plot(NA, xlim = c(0, years), ylim = c(0, 200),
       xlab = "Years", ylab = "Length (mm)",
       main = title, bty = "l")   # clean AOAS-style axes
  
  ## individual trajectories
  m <- length(res$time)
  for (i in seq_len(m)) {
    lines(res$time[[i]], res$size[[i]],
          type = "s", col = rgb(0, 0, 0, 0.3))
  }
  
  ## smoothed empirical mean
  all_times <- unlist(res$time)
  all_sizes <- unlist(res$size)
  lines(smooth.spline(all_times, all_sizes, spar = 0.7),
        col = "cyan", lwd = 3)
  
  ## theoretical VBGF curve
  L_vbgf <- Linf_pop * (1 - exp(-k * t_grid))
  lines(t_grid, L_vbgf, lwd = 2, lty = 2)
}

## 2×2 figure: fixed vs random-effects L∞
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

## 1. Females – Fixed L∞
res_f_fixed <- simulate_growth(
  n        = n,
  ni_vals  = 4:15,
  T0_range = c(10, 150),
  Linf_type = "fixed",
  Linf_par  = Linf_f_fixed,
  k = k_f, zeta = zeta_f, gam = gam_f, b1 = b1_f, b2 = b2_f
)
plot_panel(res_f_fixed, k_f, Linf_f_fixed,
           expression("Females – Fixed " * L[infinity]))

## 2. Males – Fixed L∞
res_m_fixed <- simulate_growth(
  n        = n,
  ni_vals  = 2:15,
  T0_range = c(5, 150),
  Linf_type = "fixed",
  Linf_par  = Linf_m_fixed,
  k = k_m, zeta = zeta_m, gam = gam_m, b1 = b1_m, b2 = b2_m
)
plot_panel(res_m_fixed, k_m, Linf_m_fixed,
           expression("Males – Fixed " * L[infinity]))

## 3. Females – Random-effects L∞
res_f_random <- simulate_growth(
  n        = n,
  ni_vals  = 4:15,
  T0_range = c(10, 150),
  Linf_type = "lognormal",
  Linf_par  = c(meanlog_f, sdlog_f),
  k = k_f, zeta = zeta_f, gam = gam_f, b1 = b1_f, b2 = b2_f
)
plot_panel(res_f_random, k_f, 183.27,
           expression("Females – Random-effects " * L[infinity]))

## 4. Males – Random-effects L∞
res_m_random <- simulate_growth(
  n        = n,
  ni_vals  = 2:15,
  T0_range = c(5, 150),
  Linf_type = "lognormal",
  Linf_par  = c(meanlog_m, sdlog_m),
  k = k_m, zeta = zeta_m, gam = gam_m, b1 = b1_m, b2 = b2_m
)
plot_panel(res_m_random, k_m, 184.34,
           expression("Males – Random-effects " * L[infinity]))
