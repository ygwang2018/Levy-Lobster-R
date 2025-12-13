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

##############INDEPENDENT ###########################

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



################Independent model########################
test_MIIP_ind <- function(
        Linf_f = max(subset(lobster, SEX == 1)$CL) + 10,
        Linf_m = max(subset(lobster, SEX == 2)$CL) + 30
) {
    
    library(dplyr)
    
    dat <- lobster %>%
        arrange(LOBSTER, START) %>%
        group_by(LOBSTER) %>%
        mutate(Tminus = lag(INT)) %>%
        ungroup() %>%
        filter(!is.na(Tminus), INT > 0, INC > 0)
    
    dat <- dat %>%
        group_by(SEX) %>%
        mutate(Tminus_s = as.numeric(scale(Tminus))) %>%
        ungroup()
    
    df_f <- subset(dat, SEX == 1)
    df_m <- subset(dat, SEX == 2)
    
    BL_ll <- function(par, data, Linf_fixed) {
        
        k    <- par["k"]
        zeta <- par["zeta"]
        
        eps <- 1e-8
        ll  <- 0
        
        for (i in seq_len(nrow(data))) {
            
            Lp <- data$PL[i]
            I  <- data$INC[i]
            T  <- data$INT[i] / 365
            
            mu <- 1 - exp(-k * T)
            mu <- min(max(mu, eps), 1 - eps)
            
            x <- I / (Linf_fixed - Lp)
            x <- min(max(x, eps), 1 - eps)
            
            a <- (zeta - 1) * mu
            b <- (zeta - 1) * (1 - mu)
            
            if (a <= 0 || b <= 0) return(1e10)
            
            ll <- ll + dbeta(x, a, b, log = TRUE)
        }
        
        -ll
    }
    
    fit_BL_f <- optim(
        c(k = 0.3, zeta = 50), BL_ll,
        data = df_f, Linf_fixed = Linf_f,
        method = "L-BFGS-B",
        lower = c(1e-4, 2),
        upper = c(0.5, 300)
    )
    
    fit_BL_m <- optim(
        c(k = 0.3, zeta = 50), BL_ll,
        data = df_m, Linf_fixed = Linf_m,
        method = "L-BFGS-B",
        lower = c(1e-4, 2),
        upper = c(0.5, 300)
    )
    
    ref_f <- fit_BL_f$par
    ref_m <- fit_BL_m$par
    
        MIIP_ll <- function(par, data, Linf_fixed) {
        
        k     <- par["k"]
        zeta  <- par["zeta"]
        phi   <- par["phi"]
        beta0 <- par["beta0"]
        beta1 <- par["beta1"]
        beta2 <- par["beta2"]
        
        eps <- 1e-8
        ll  <- 0
        
        for (i in seq_len(nrow(data))) {
            
            Lp <- data$PL[i]
            I  <- data$INC[i]
            T  <- data$INT[i] / 365
            
            ## ---- MI: BL ----
            mu <- 1 - exp(-k * T)
            mu <- min(max(mu, eps), 1 - eps)
            
            x <- I / (Linf_fixed - Lp)
            x <- min(max(x, eps), 1 - eps)
            
            a <- (zeta - 1) * mu
            b <- (zeta - 1) * (1 - mu)
            
            if (a <= 0 || b <= 0) return(1e10)
            
            ll <- ll + dbeta(x, a, b, log = TRUE)
            
            ## ---- IP: Gamma ----
            mu_T <- exp(beta0 +
                            beta1 * Lp +
                            beta2 * data$Tminus_s[i])
            
            if (mu_T <= 0) return(1e10)
            
            ll <- ll + dgamma(
                data$INT[i],
                shape = phi,
                scale = mu_T / phi,
                log = TRUE
            )
        }
        
        -ll
    }
    
    init  <- c(k=0.3, zeta=50, phi=5, beta0=1, beta1=0.01, beta2=0.2)
    lower <- c(1e-4, 2, 0.5, -10, -0.05, -5)
    upper <- c(0.5, 300, 50, 10, 0.05, 5)
    
    fit_f <- optim(
        init, MIIP_ll,
        data = df_f, Linf_fixed = Linf_f,
        method = "L-BFGS-B",
        lower = lower, upper = upper,
        hessian = TRUE
    )
    
    fit_m <- optim(
        init, MIIP_ll,
        data = df_m, Linf_fixed = Linf_m,
        method = "L-BFGS-B",
        lower = lower, upper = upper,
        hessian = TRUE
    )
    
    SE_f <- sqrt(mean(diag(solve(fit_f$hessian)), na.rm = TRUE))
    SE_m <- sqrt(mean(diag(solve(fit_m$hessian)), na.rm = TRUE))
    
    Bias_f <- sqrt(mean((fit_f$par[c("k","zeta")] - ref_f)^2))
    Bias_m <- sqrt(mean((fit_m$par[c("k","zeta")] - ref_m)^2))
    
    RMSE_f <- sqrt(SE_f^2 + Bias_f^2)
    RMSE_m <- sqrt(SE_m^2 + Bias_m^2)
    
    results <- data.frame(
        Sex   = c("Female","Male"),
        Linf  = c(Linf_f, Linf_m),
        k     = c(fit_f$par["k"], fit_m$par["k"]),
        SE    = c(SE_f, SE_m),
        Bias  = c(Bias_f, Bias_m),
        RMSE  = c(RMSE_f, RMSE_m),
        NegLL = c(fit_f$value, fit_m$value),
        Converged = c(
            fit_f$convergence == 0,
            fit_m$convergence == 0
        )
    )
    
    print(results)
    invisible(results)
}
