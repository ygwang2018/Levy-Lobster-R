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

laplace_logint <- function(f, lo, hi) {

  opt <- try(optimize(function(u) -f(u), c(lo, hi)), silent = TRUE)
  if (inherits(opt, "try-error")) return(-Inf)

  u0 <- opt$minimum
  f0 <- f(u0)
  if (!is.finite(f0)) return(-Inf)

  h <- 1e-4
  fpp <- (f(u0 + h) - 2 * f0 + f(u0 - h)) / h^2

  if (!is.finite(fpp) || fpp >= 0) return(-Inf)

  f0 + 0.5 * log(2 * pi) - 0.5 * log(-fpp)
}

simulate_ind_MI <- function(
  n_ind = 80,
  n_obs = 10,
  true
) {

  dat <- vector("list", n_ind * n_obs)
  idx <- 1
  id  <- 1

  for (i in seq_len(n_ind)) {

    Linf_i <- rlnorm(
      1,
      meanlog = log(true$Linf_mean) - 0.5 * true$Linf_sdlog^2,
      sdlog   = true$Linf_sdlog
    )

    PL     <- runif(1, 40, 80)
    Tminus <- runif(1, 0.5, 2)

    for (j in seq_len(n_obs)) {

      ## FIXED INT DESIGN
      INT <- runif(1, 0.5, 2.5)

      ## MI: Beta–Lognormal
      mu <- 1 - exp(-true$k * INT)
      mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)

      a <- (true$zeta - 1) * mu
      b <- (true$zeta - 1) * (1 - mu)

      x   <- rbeta(1, a, b)
      INC <- x * (Linf_i - PL)

      dat[[idx]] <- data.frame(
        LOBSTER = id,
        PL      = PL,
        INC     = INC,
        INT     = INT
      )
      idx <- idx + 1

      PL     <- PL + INC
      Tminus <- INT
    }

    id <- id + 1
  }

  do.call(rbind, dat)
}


fit_ind_MI <- function(dat) {

  nll <- function(par) {

    k     <- par["k"]
    zeta  <- par["zeta"]
    mu_L  <- par["mu_L"]
    sd_L  <- par["sd_L"]

    if (k <= 0 || zeta <= 1 || sd_L <= 0)
      return(1e10)

    eps <- 1e-8
    ll  <- 0

    for (id in unique(dat$LOBSTER)) {

      d <- dat[dat$LOBSTER == id, ]
      maxPL <- max(d$PL)

      f_u <- function(u) {

        Linf <- exp(u)
        if (Linf <= maxPL) return(-Inf)

        ll_i <- dnorm(u, mu_L, sd_L, log = TRUE)

        for (j in seq_len(nrow(d))) {

          mu <- 1 - exp(-k * d$INT[j])
          mu <- pmin(pmax(mu, eps), 1 - eps)

          denom <- Linf - d$PL[j]
          if (denom <= eps) return(-Inf)

          x <- d$INC[j] / denom
          if (x <= eps || x >= 1 - eps) return(-Inf)

          a <- (zeta - 1) * mu
          b <- (zeta - 1) * (1 - mu)

          ll_i <- ll_i + dbeta(x, a, b, log = TRUE)
        }

        ll_i
      }

      lo <- log(maxPL + 1)
      hi <- log(maxPL + 3 * exp(mu_L + 2 * sd_L))

      logLi <- laplace_logint(f_u, lo, hi)
      if (!is.finite(logLi)) return(1e10)

      ll <- ll + logLi
    }

    -ll
  }

  init <- c(
    k     = 0.25,
    zeta  = 120,
    mu_L  = log(mean(dat$PL) + 80),
    sd_L  = 0.08
  )

  lower <- c(0.05, 50, log(60), 0.02)
  upper <- c(0.8, 400, log(400), 0.4)

  optim(init, nll, method = "L-BFGS-B",
        lower = lower, upper = upper,
        control = list(maxit = 1500))
}


extract_k_Linf <- function(fit) {
  c(
    Linf = exp(fit$par["mu_L"] + 0.5 * fit$par["sd_L"]^2),
    k    = fit$par["k"]
  )
}

## Monte Carlo (sex-specific)
R <- 1000

true_F <- list(
  Linf_mean = 180,
  Linf_sdlog = 0.08,
  k = 0.28,
  zeta = 120
)

true_M <- list(
  Linf_mean = 210,
  Linf_sdlog = 0.08,
  k = 0.24,
  zeta = 120
)

storeF <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("Linf","k")))
storeM <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("Linf","k")))

for (r in seq_len(R)) {

  datF <- simulate_ind_MI(true = true_F)
  datM <- simulate_ind_MI(true = true_M)

  fitF <- try(fit_ind_MI(datF), silent = TRUE)
  fitM <- try(fit_ind_MI(datM), silent = TRUE)

  if (!inherits(fitF, "try-error") && fitF$convergence == 0)
    storeF[r, ] <- extract_k_Linf(fitF)

  if (!inherits(fitM, "try-error") && fitM$convergence == 0)
    storeM[r, ] <- extract_k_Linf(fitM)

  if (r %% 50 == 0)
    cat("Completed", r, "\n")
}

## Monte Carlo summary

MC_summary <- function(est, truth) {

  est <- est[complete.cases(est), , drop = FALSE]

  data.frame(
    Parameter = colnames(est),
    Mean  = colMeans(est),
    SE    = apply(est, 2, sd),
    Bias  = colMeans(est) - truth[colnames(est)],
    RMSE  = sqrt(colMeans((est - truth[colnames(est)])^2))
  )
}

print(MC_summary(storeF, c(Linf = 180, k = 0.28)))
print(MC_summary(storeM, c(Linf = 210, k = 0.24)))
