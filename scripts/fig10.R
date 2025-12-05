ln_params <- function(mean, sd) {
  var <- sd^2
  meanlog <- log(mean^2 / sqrt(var + mean^2))
  sdlog   <- sqrt(log(1 + var / mean^2))
  return(c(meanlog = meanlog, sdlog = sdlog))
}

simulate_growth <- function(
    n, ni_vals, T0_range,
    Linf_type = "fixed", Linf_par,
    k, zeta, gam, b1, b2
) {
  res <- vector("list", n)
  
  for (i in 1:n) {
    ni <- sample(ni_vals, 1)
    T0 <- runif(1, T0_range[1], T0_range[2])
    Ti <- cumsum(rexp(ni, rate = 1/365))
    
    # Generate Linf
    Linf <- if (Linf_type == "fixed") {
      Linf_par
    } else if (Linf_type == "lognormal") {
      rlnorm(1, meanlog = Linf_par[1], sdlog = Linf_par[2])
    }
    
    # Expected increments
    E_MI <- (Linf - b2) * (1 - exp(-k * (Ti - T0)))
    
    # Gamma variability
    shape <- zeta
    rate  <- zeta / E_MI
    MI    <- rgamma(ni, shape = shape, rate = rate)
    
    L <- b2 + cumsum(MI)
    
    res[[i]] <- data.frame(
      i = i, ni = ni, Ti = Ti, L = L, Linf = Linf
    )
  }
  
  return(do.call(rbind, res))
}

plot_panel <- function(df, k, Linf_mean, title_expr) {
  plot(df$Ti/365, df$L,
       type = "n",
       xlab = "Years", ylab = "Carapace Length (mm)",
       main = title_expr)
  
  ids <- unique(df$i)
  for (id in ids) {
    d <- df[df$i == id, ]
    lines(d$Ti/365, d$L, col = rgb(0,0,0,0.25))
  }
  
  # Expected VB curve
  t_grid <- seq(0, 10, length.out = 300)
  Lhat <- Linf_mean * (1 - exp(-k * t_grid))
  lines(t_grid, Lhat, col = "red", lwd = 2)
}

test_fig10 <- function() {
  
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
  
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  png("results/figures/growth_panel.png", width = 1200, height = 900, res = 150)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # 1. Females – Fixed L∞
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
  
  # 2. Males – Fixed L∞
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
  
  # 3. Females – Random-effects L∞
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
  
  # 4. Males – Random-effects L∞
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
  
  dev.off()
  
  cat("growth_panel.png saved to results/figures.\n")
  
  invisible(list(fixed_f = res_f_fixed,
                 fixed_m = res_m_fixed,
                 random_f = res_f_random,
                 random_m = res_m_random))
}
