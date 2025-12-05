library(dplyr)

# 1. Convert (mean, sd) to (meanlog, sdlog)
ln_params <- function(mean, sd) {
  var <- sd^2
  meanlog <- log(mean^2 / sqrt(var + mean^2))
  sdlog   <- sqrt(log(1 + var / mean^2))
  return(c(meanlog = meanlog, sdlog = sdlog))
}

# 2. Correct von Bertalanffy–Gamma growth simulator
simulate_growth <- function(
    n, ni_vals, T0_range,
    Linf_type = "fixed", Linf_par,
    k, L0 = 5           # starting length = 5 mm
) {
  
  out <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    ni  <- sample(ni_vals, 1)
    
    # Exponential inter-moult times
    dT <- rexp(ni, rate = 1/365)
    Ti <- cumsum(dT)
    
    # Generate individual asymptotic length
    Linf <- if (Linf_type == "fixed") {
      Linf_par
    } else {
      rlnorm(1, meanlog = Linf_par[1], sdlog = Linf_par[2])
    }
    
    # Track actual growth
    L <- numeric(ni)
    L_prev <- L0
    
    for (j in seq_len(ni)) {
      # Expected VB increment
      mu <- (Linf - L_prev) * (1 - exp(-k * dT[j]))
      
      # Gamma noise with moderate shape → realistic variability
      shape <- 15
      scale <- mu / shape
      MI <- rgamma(1, shape = shape, scale = scale)
      
      L[j] <- L_prev + MI
      L_prev <- L[j]
    }
    
    out[[i]] <- data.frame(
      id   = i,
      Ti   = Ti,
      L    = L,
      Linf = Linf
    )
  }
  
  bind_rows(out)
}

# 3. Plot panel (fixed or random-effects L∞)
plot_panel <- function(df, k, Linf_mean, title_expr, L0 = 5) {
  
  plot(df$Ti/365, df$L,
       type="n",
       xlab="Years", ylab="Length (mm)",
       main=title_expr)
  
  # Individual growth trajectories
  ids <- unique(df$id)
  for (id in ids) {
    d <- df[df$id == id, ]
    lines(d$Ti/365, d$L, col = rgb(0,0,0,0.25), lwd=1)
  }
  
  # Correct VB curve
  t <- seq(0, max(df$Ti)/365, length.out=300)
  Lhat <- Linf_mean - (Linf_mean - L0)*exp(-k*t)
  lines(t, Lhat, col="cyan", lwd=2)
}

# 4. Main figure wrapper
test_fig10 <- function() {
  
  set.seed(123)
  
  n <- 80
  
  # === Female parameters ===
  Linf_f_fixed <- 183.87
  f_ln <- ln_params(183.27, 12)
  k_f  <- 0.283
  
  # === Male parameters ===
  Linf_m_fixed <- 228.04
  m_ln <- ln_params(184.34, 12)
  k_m  <- 0.2473
  
    
  png("results/figures/growth_panel.png",
      width = 1400, height = 1100, res = 150)
  
  par(mfrow = c(2, 2), mar = c(4,4,3,1))
  
  # 1. Females – Fixed L∞
  res_f_fixed <- simulate_growth(
    n = n, ni_vals = 4:15,
    T0_range = c(10,150),
    Linf_type = "fixed",
    Linf_par  = Linf_f_fixed,
    k = k_f
  )
  plot_panel(res_f_fixed, k_f, Linf_f_fixed,
             expression("Females – Fixed " * L[infinity]))
  
  # 2. Males – Fixed L∞
  res_m_fixed <- simulate_growth(
    n = n, ni_vals = 2:15,
    T0_range = c(5,150),
    Linf_type = "fixed",
    Linf_par  = Linf_m_fixed,
    k = k_m
  )
  plot_panel(res_m_fixed, k_m, Linf_m_fixed,
             expression("Males – Fixed " * L[infinity]))
  
  # 3. Females – Random-effects L∞
  res_f_random <- simulate_growth(
    n = n, ni_vals = 4:15,
    T0_range = c(10,150),
    Linf_type = "lognormal",
    Linf_par  = c(f_ln["meanlog"], f_ln["sdlog"]),
    k = k_f
  )
  plot_panel(res_f_random, k_f, 183.27,
             expression("Females – Random-effects " * L[infinity]))
  
  # 4. Males – Random-effects L∞
  res_m_random <- simulate_growth(
    n = n, ni_vals = 2:15,
    T0_range = c(5,150),
    Linf_type = "lognormal",
    Linf_par  = c(m_ln["meanlog"], m_ln["sdlog"]),
    k = k_m
  )
  plot_panel(res_m_random, k_m, 184.34,
             expression("Males – Random-effects " * L[infinity]))
  
  dev.off()
  }

