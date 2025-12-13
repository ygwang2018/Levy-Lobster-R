##########Lognormal Linf###########################
simulate_lognormal <- function(R = 1000, n_f = 100, n_m = 100) {
  
  set.seed(123)
  
  k_range_f    <- c(0.20, 0.50)
  Linf_range_f <- c(150, 230)
  
  k_range_m    <- c(0.2, 0.5)
  Linf_range_m <- c(150, 230)
  
  sigma_eps <- 1.5
  sdlog_L   <- 0.05
  
  ## Observation times
  t_obs <- sort(runif(6, 0.5, 10))
  
  sim_data <- function(n, k_true, Linf_mean) {
    
    mu_logn <- log(Linf_mean) - 0.5 * sdlog_L^2
    Linf_i  <- rlnorm(n, mu_logn, sdlog_L)
    
    do.call(
      rbind,
      lapply(seq_len(n), function(i) {
        data.frame(
          id = i,
          t  = t_obs,
          L  = Linf_i[i] * (1 - exp(-k_true * t_obs)) +
            rnorm(length(t_obs), 0, sigma_eps)
        )
      })
    )
  }
  
  
  nll <- function(par, dat) {
    
    k    <- par[1]
    Linf <- par[2]
    
    if (k <= 0 || Linf <= max(dat$L)) return(1e10)
    
    mu <- Linf * (1 - exp(-k * dat$t))
    -sum(dnorm(dat$L, mu, sigma_eps, log = TRUE))
  }
  
  khat_f <- Linfhat_f <- numeric(R)
  khat_m <- Linfhat_m <- numeric(R)
  ktrue_f <- Ltrue_f <- numeric(R)
  ktrue_m <- Ltrue_m <- numeric(R)
  
  
  for (r in seq_len(R)) {
    
    ## --- FEMALE ---
    kf <- runif(1, k_range_f[1], k_range_f[2])
    Lf <- runif(1, Linf_range_f[1], Linf_range_f[2])
    
    dat_f <- sim_data(n_f, kf, Lf)
    
    opt_f <- optim(
      par    = c(mean(k_range_f), mean(Linf_range_f)),
      fn     = nll,
      dat    = dat_f,
      method = "L-BFGS-B",
      lower  = c(k_range_f[1], Linf_range_f[1]),
      upper  = c(k_range_f[2], Linf_range_f[2])
    )
    
    khat_f[r]    <- opt_f$par[1]
    Linfhat_f[r] <- opt_f$par[2]
    ktrue_f[r]   <- kf
    Ltrue_f[r]   <- Lf
    
    ## --- MALE ---
    km <- runif(1, k_range_m[1], k_range_m[2])
    Lm <- runif(1, Linf_range_m[1], Linf_range_m[2])
    
    dat_m <- sim_data(n_m, km, Lm)
    
    opt_m <- optim(
      par    = c(mean(k_range_m), mean(Linf_range_m)),
      fn     = nll,
      dat    = dat_m,
      method = "L-BFGS-B",
      lower  = c(k_range_m[1], Linf_range_m[1]),
      upper  = c(k_range_m[2], Linf_range_m[2])
    )
    
    khat_m[r]    <- opt_m$par[1]
    Linfhat_m[r] <- opt_m$par[2]
    ktrue_m[r]   <- km
    Ltrue_m[r]   <- Lm
  }
  
  summary_sex <- rbind(
    data.frame(
      Sex = "Female",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_f), mean(Linfhat_f)),
      Bias = c(mean(khat_f - ktrue_f), mean(Linfhat_f - Ltrue_f)),
      SE   = c(sd(khat_f), sd(Linfhat_f)),
      RMSE = c(
        sqrt(mean((khat_f - ktrue_f)^2)),
        sqrt(mean((Linfhat_f - Ltrue_f)^2))
      )
    ),
    data.frame(
      Sex = "Male",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_m), mean(Linfhat_m)),
      Bias = c(mean(khat_m - ktrue_m), mean(Linfhat_m - Ltrue_m)),
      SE   = c(sd(khat_m), sd(Linfhat_m)),
      RMSE = c(
        sqrt(mean((khat_m - ktrue_m)^2)),
        sqrt(mean((Linfhat_m - Ltrue_m)^2))
      )
    )
  )
  
  print(summary_sex)
}
simulate_lognormal()

############Fixed Linf#############################
simulate_fixed <- function(R = 1000, n_f = 100, n_m = 100) {
  
  set.seed(123)
  
  k_range_f    <- c(0.20, 0.45)
  Linf_range_f <- c(150, 250)
  
  k_range_m    <- c(0.2, 0.50)
  Linf_range_m <- c(150, 250)
  
  sigma_eps <- 1.5
  
  ## Observation times
  t_obs <- sort(runif(6, 0.5, 10))
  
    sim_data <- function(n, k_true, Linf_true) {
    
    do.call(
      rbind,
      lapply(seq_len(n), function(i) {
        data.frame(
          id = i,
          t  = t_obs,
          L  = Linf_true * (1 - exp(-k_true * t_obs)) +
            rnorm(length(t_obs), 0, sigma_eps)
        )
      })
    )
  }
  
  nll <- function(par, dat) {
    
    k    <- par[1]
    Linf <- par[2]
    
    if (k <= 0 || Linf <= max(dat$L)) return(1e10)
    
    mu <- Linf * (1 - exp(-k * dat$t))
    -sum(dnorm(dat$L, mu, sigma_eps, log = TRUE))
  }
  
  khat_f <- Linfhat_f <- numeric(R)
  khat_m <- Linfhat_m <- numeric(R)
  
  ktrue_f <- Ltrue_f <- numeric(R)
  ktrue_m <- Ltrue_m <- numeric(R)
  
    for (r in seq_len(R)) {
    
    ## ---------- FEMALE ----------
    kf <- runif(1, k_range_f[1], k_range_f[2])
    Lf <- runif(1, Linf_range_f[1], Linf_range_f[2])
    
    dat_f <- sim_data(n_f, kf, Lf)
    
    opt_f <- optim(
      par    = c(mean(k_range_f), mean(Linf_range_f)),
      fn     = nll,
      dat    = dat_f,
      method = "L-BFGS-B",
      lower  = c(k_range_f[1], Linf_range_f[1]),
      upper  = c(k_range_f[2], Linf_range_f[2])
    )
    
    khat_f[r]    <- opt_f$par[1]
    Linfhat_f[r] <- opt_f$par[2]
    ktrue_f[r]   <- kf
    Ltrue_f[r]   <- Lf
    
    ## ---------- MALE ----------
    km <- runif(1, k_range_m[1], k_range_m[2])
    Lm <- runif(1, Linf_range_m[1], Linf_range_m[2])
    
    dat_m <- sim_data(n_m, km, Lm)
    
    opt_m <- optim(
      par    = c(mean(k_range_m), mean(Linf_range_m)),
      fn     = nll,
      dat    = dat_m,
      method = "L-BFGS-B",
      lower  = c(k_range_m[1], Linf_range_m[1]),
      upper  = c(k_range_m[2], Linf_range_m[2])
    )
    
    khat_m[r]    <- opt_m$par[1]
    Linfhat_m[r] <- opt_m$par[2]
    ktrue_m[r]   <- km
    Ltrue_m[r]   <- Lm
  }
  
    summary_fixed_sex <- rbind(
    data.frame(
      Sex = "Female",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_f), mean(Linfhat_f)),
      Bias = c(mean(khat_f - ktrue_f),
               mean(Linfhat_f - Ltrue_f)),
      SE   = c(sd(khat_f), sd(Linfhat_f)),
      RMSE = c(
        sqrt(mean((khat_f - ktrue_f)^2)),
        sqrt(mean((Linfhat_f - Ltrue_f)^2))
      )
    ),
    data.frame(
      Sex = "Male",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_m), mean(Linfhat_m)),
      Bias = c(mean(khat_m - ktrue_m),
               mean(Linfhat_m - Ltrue_m)),
      SE   = c(sd(khat_m), sd(Linfhat_m)),
      RMSE = c(
        sqrt(mean((khat_m - ktrue_m)^2)),
        sqrt(mean((Linfhat_m - Ltrue_m)^2))
      )
    )
  )
  
  print(summary_fixed_sex)
}

simulate_fixed()

########################Gamma Linf##############################
simulate_gamma<- function(R = 1000, n_f = 100, n_m = 100) {
  
  set.seed(123)
  
  k_range_f    <- c(0.20, 0.5)
  Linf_range_f <- c(150, 250)
  
  k_range_m    <- c(0.2, 0.50)
  Linf_range_m <- c(150, 250)
  
  sigma_eps <- 1.5
  
  ## Gamma dispersion (same for both sexes)
  shape_L <- 200   # larger = less dispersion
  
  ## Observation times
  t_obs <- sort(runif(6, 0.5, 10))
  
  sim_data <- function(n, k_true, Linf_mean) {
    
    scale_L <- Linf_mean / shape_L
    Linf_i  <- rgamma(n, shape = shape_L, scale = scale_L)
    
    do.call(
      rbind,
      lapply(seq_len(n), function(i) {
        data.frame(
          id = i,
          t  = t_obs,
          L  = Linf_i[i] * (1 - exp(-k_true * t_obs)) +
            rnorm(length(t_obs), 0, sigma_eps)
        )
      })
    )
  }
  
  nll <- function(par, dat) {
    
    k    <- par[1]
    Linf <- par[2]
    
    if (k <= 0 || Linf <= max(dat$L)) return(1e10)
    
    mu <- Linf * (1 - exp(-k * dat$t))
    -sum(dnorm(dat$L, mu, sigma_eps, log = TRUE))
  }
  
  khat_f <- Linfhat_f <- numeric(R)
  khat_m <- Linfhat_m <- numeric(R)
  
  ktrue_f <- Ltrue_f <- numeric(R)
  ktrue_m <- Ltrue_m <- numeric(R)
  
  for (r in seq_len(R)) {
    
    ## ---------- FEMALE ----------
    kf <- runif(1, k_range_f[1], k_range_f[2])
    Lf <- runif(1, Linf_range_f[1], Linf_range_f[2])
    
    dat_f <- sim_data(n_f, kf, Lf)
    
    opt_f <- optim(
      par    = c(mean(k_range_f), mean(Linf_range_f)),
      fn     = nll,
      dat    = dat_f,
      method = "L-BFGS-B",
      lower  = c(k_range_f[1], Linf_range_f[1]),
      upper  = c(k_range_f[2], Linf_range_f[2])
    )
    
    khat_f[r]    <- opt_f$par[1]
    Linfhat_f[r] <- opt_f$par[2]
    ktrue_f[r]   <- kf
    Ltrue_f[r]   <- Lf
    
    ## ---------- MALE ----------
    km <- runif(1, k_range_m[1], k_range_m[2])
    Lm <- runif(1, Linf_range_m[1], Linf_range_m[2])
    
    dat_m <- sim_data(n_m, km, Lm)
    
    opt_m <- optim(
      par    = c(mean(k_range_m), mean(Linf_range_m)),
      fn     = nll,
      dat    = dat_m,
      method = "L-BFGS-B",
      lower  = c(k_range_m[1], Linf_range_m[1]),
      upper  = c(k_range_m[2], Linf_range_m[2])
    )
    
    khat_m[r]    <- opt_m$par[1]
    Linfhat_m[r] <- opt_m$par[2]
    ktrue_m[r]   <- km
    Ltrue_m[r]   <- Lm
  }
  
  summary_gamma_sex <- rbind(
    data.frame(
      Sex = "Female",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_f), mean(Linfhat_f)),
      Bias = c(mean(khat_f - ktrue_f),
               mean(Linfhat_f - Ltrue_f)),
      SE   = c(sd(khat_f), sd(Linfhat_f)),
      RMSE = c(
        sqrt(mean((khat_f - ktrue_f)^2)),
        sqrt(mean((Linfhat_f - Ltrue_f)^2))
      )
    ),
    data.frame(
      Sex = "Male",
      Parameter = c("k", "Linf"),
      Mean = c(mean(khat_m), mean(Linfhat_m)),
      Bias = c(mean(khat_m - ktrue_m),
               mean(Linfhat_m - Ltrue_m)),
      SE   = c(sd(khat_m), sd(Linfhat_m)),
      RMSE = c(
        sqrt(mean((khat_m - ktrue_m)^2)),
        sqrt(mean((Linfhat_m - Ltrue_m)^2))
      )
    )
  )
  
  print(summary_gamma_sex)
}
simulate_gamma()
