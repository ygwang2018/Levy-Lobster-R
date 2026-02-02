set.seed(123)

R_mc <- 1000

n_f <- 100
n_m <- 100

true_f <- list(k = 0.28, Linf = 185.32)
true_m <- list(k = 0.30, Linf = 201.48)

t_obs_fixed <- c(0.4, 2.0)                         # 2-point design
t_obs_long  <- c(0.3, 0.8, 1.5, 2.5, 3.0, 4.0, 4.5) # 7-point design

sigma_eps <- 11.5

sdlog_L <- 0.08     # Lognormal sd(log L_infty)
shape_L <- 130      # Gamma shape (CV ≈ 1/sqrt(shape) ≈ 0.088)

lower_bounds <- c(0.05, 50)   # (k, Linf)
upper_bounds <- c(1.00, 400)
opt_control  <- list(maxit = 5000)

tol_bound <- 1e-6   # detect exact boundary hits
tol_clip  <- 1e-8   # numerical eps

make_nll_vbgf <- function(sigma_eps) {
  function(par, dat) {
    k    <- par[1]
    Linf <- par[2]
    # hard constraints
    if (k <= 0 || Linf <= max(dat$L)) return(1e10)
    mu <- Linf * (1 - base::exp(-k * dat$t))
    # negative log-likelihood
    -base::sum(stats::dnorm(dat$L, mean = mu, sd = sigma_eps, log = TRUE))
  }
}

.clamp_open <- function(x, lo, hi, eps = 1e-6) {
  x <- pmax(x, lo + eps)
  x <- pmin(x, hi - eps)
  x
}

fit_one_multistart <- function(dat, k_init, Linf_init, sigma_eps, n_starts = 5) {
  nll <- make_nll_vbgf(sigma_eps)
  
  starts <- matrix(NA_real_, n_starts, 2)
  starts[1, ] <- c(k_init, Linf_init)
  # jittered starts
  for (s in 2:n_starts) {
    k0 <- k_init  * (1 + stats::rnorm(1, 0, 0.25))
    L0 <- Linf_init * (1 + stats::rnorm(1, 0, 0.20))
    starts[s, ] <- c(
      .clamp_open(k0, lower_bounds[1], upper_bounds[1]),
      .clamp_open(L0, lower_bounds[2], upper_bounds[2])
    )
  }
  
  best <- NULL
  for (s in 1:n_starts) {
    res <- try(
      stats::optim(
        par    = starts[s, ],
        fn     = nll,
        dat    = dat,
        method = "L-BFGS-B",
        lower  = lower_bounds,
        upper  = upper_bounds,
        control= opt_control
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error")) next
    if (is.null(best) || res$value < best$value) best <- res
  }
  best
}

sim_fixed_dat <- function(N, k, L, t_obs, sigma_eps) {
  nt    <- base::length(t_obs)
  t_all <- base::rep(t_obs, N)
  id    <- base::rep(1:N, each = nt)
  Lobs  <- L * (1 - base::exp(-k * t_all)) + stats::rnorm(N*nt, 0, sigma_eps)
  base::data.frame(id = id, t = t_all, L = Lobs)
}

sim_logn_dat <- function(N, k, Linf_mean, t_obs, sigma_eps, sdlog_L) {
  nt      <- base::length(t_obs)
  mu_logn <- base::log(Linf_mean) - 0.5 * sdlog_L^2
  Linf_i  <- stats::rlnorm(N, mu_logn, sdlog_L)
  t_all   <- base::rep(t_obs, N)
  id      <- base::rep(1:N, each = nt)
  Linf_rep<- base::rep(Linf_i, each = nt)
  Lobs    <- Linf_rep * (1 - base::exp(-k*t_all)) + stats::rnorm(N*nt, 0, sigma_eps)
  base::data.frame(id = id, t = t_all, L = Lobs)
}

sim_gamma_dat <- function(N, k, Linf_mean, t_obs, sigma_eps, shape_L) {
  nt      <- base::length(t_obs)
  Linf_i  <- stats::rgamma(N, shape = shape_L, scale = Linf_mean/shape_L)
  t_all   <- base::rep(t_obs, N)
  id      <- base::rep(1:N, each = nt)
  Linf_rep<- base::rep(Linf_i, each = nt)
  Lobs    <- Linf_rep * (1 - base::exp(-k*t_all)) + stats::rnorm(N*nt, 0, sigma_eps)
  base::data.frame(id = id, t = t_all, L = Lobs)
}


sim_fun_fixed <- function(N, k, Linf, t_obs, sigma_eps, sdlog_L = NULL, shape_L = NULL) {
  sim_fixed_dat(N, k, Linf, t_obs, sigma_eps)
}

sim_fun_logn <- function(N, k, Linf, t_obs, sigma_eps, sdlog_L = sdlog_L, shape_L = NULL) {
  sim_logn_dat(N, k, Linf, t_obs, sigma_eps, sdlog_L)
}

sim_fun_gamma <- function(N, k, Linf, t_obs, sigma_eps, sdlog_L = NULL, shape_L = shape_L) {
  sim_gamma_dat(N, k, Linf, t_obs, sigma_eps, shape_L)
}

summarise_vec <- function(x) {
  if (length(x) < 2) {
    c(mean = mean(x), sd = 0, q025 = x, q975 = x, n = length(x))
  } else {
    c(mean = mean(x), sd = stats::sd(x),
      q025 = stats::quantile(x, 0.025),
      q975 = stats::quantile(x, 0.975),
      n = length(x))
  }
}

mc_summary <- function(khat, Lhat, k_true, L_true,
                       k_bound, L_bound, conv_code) {
  
  keep <- which(conv_code == 0 & !k_bound & !L_bound)
  kh <- khat[keep]; Lh <- Lhat[keep]
  
  if (length(kh) < 5) {
    keep <- which(conv_code == 0 & !k_bound)
    kh <- khat[keep]; Lh <- Lhat[keep]
  }
  if (length(kh) < 5) {
    keep <- which(conv_code == 0)
    kh <- khat[keep]; Lh <- Lhat[keep]
  }
  
  k_stats <- summarise_vec(kh)
  L_stats <- summarise_vec(Lh)
  
  base::data.frame(
    Parameter = c("k","Linf"),
    True      = c(k_true, L_true),
    Mean      = c(k_stats["mean"], L_stats["mean"]),
    Bias      = c(k_stats["mean"] - k_true, L_stats["mean"] - L_true),
    SE        = c(k_stats["sd"],   L_stats["sd"]),
    RMSE      = c(base::sqrt(base::mean((kh - k_true)^2)),
                  base::sqrt(base::mean((Lh - L_true)^2))),
    CI_low    = c(k_stats["q025"], L_stats["q025"]),
    CI_high   = c(k_stats["q975"], L_stats["q975"]),
    N_used    = c(k_stats["n"],    L_stats["n"])
  )
}

run_mc <- function(R, sim_fun, t_obs, sigma_eps,
                   sdlog_L = sdlog_L, shape_L = shape_L,
                   trueF = true_f, trueM = true_m, label = "") {
  
  kF <- LF <- base::numeric(R)
  kM <- LM <- base::numeric(R)
  
  convF <- convM <- base::integer(R)
  kB_F <- LB_F <- kB_M <- LB_M <- base::logical(R)
  
  for (r in 1:R) {
    # Female
    datF <- sim_fun(
      N = n_f, k = trueF$k, Linf = trueF$Linf,
      t_obs = t_obs, sigma_eps = sigma_eps,
      sdlog_L = sdlog_L, shape_L = shape_L
    )
    fitF <- fit_one_multistart(datF, trueF$k, trueF$Linf, sigma_eps, n_starts = 5)
    if (!is.null(fitF)) {
      kF[r]   <- fitF$par[1]
      LF[r]   <- fitF$par[2]
      convF[r]<- fitF$convergence
      kB_F[r] <- abs(kF[r] - lower_bounds[1]) < tol_bound |
        abs(kF[r] - upper_bounds[1]) < tol_bound
      LB_F[r] <- abs(LF[r] - lower_bounds[2]) < tol_bound |
        abs(LF[r] - upper_bounds[2]) < tol_bound
    } else {
      kF[r] <- NA_real_; LF[r] <- NA_real_; convF[r] <- 99L; kB_F[r] <- TRUE; LB_F[r] <- TRUE
    }
    
    # Male
    datM <- sim_fun(
      N = n_m, k = trueM$k, Linf = trueM$Linf,
      t_obs = t_obs, sigma_eps = sigma_eps,
      sdlog_L = sdlog_L, shape_L = shape_L
    )
    fitM <- fit_one_multistart(datM, trueM$k, trueM$Linf, sigma_eps, n_starts = 5)
    if (!is.null(fitM)) {
      kM[r]   <- fitM$par[1]
      LM[r]   <- fitM$par[2]
      convM[r]<- fitM$convergence
      kB_M[r] <- abs(kM[r] - lower_bounds[1]) < tol_bound |
        abs(kM[r] - upper_bounds[1]) < tol_bound
      LB_M[r] <- abs(LM[r] - lower_bounds[2]) < tol_bound |
        abs(LM[r] - upper_bounds[2]) < tol_bound
    } else {
      kM[r] <- NA_real_; LM[r] <- NA_real_; convM[r] <- 99L; kB_M[r] <- TRUE; LB_M[r] <- TRUE
    }
    
    if (r %% 100 == 0) base::message(label, " : ", r, "/", R)
  }
  
  # Diagnostics
  diag_f <- base::data.frame(
    Sex            = "Female",
    Scenario       = label,
    Converged      = base::sum(convF == 0, na.rm = TRUE),
    Total          = R,
    k_boundary     = base::sum(kB_F, na.rm = TRUE),
    Linf_boundary  = base::sum(LB_F, na.rm = TRUE),
    NA_fits        = base::sum(!base::is.finite(kF) | !base::is.finite(LF))
  )
  diag_m <- base::data.frame(
    Sex            = "Male",
    Scenario       = label,
    Converged      = base::sum(convM == 0, na.rm = TRUE),
    Total          = R,
    k_boundary     = base::sum(kB_M, na.rm = TRUE),
    Linf_boundary  = base::sum(LB_M, na.rm = TRUE),
    NA_fits        = base::sum(!base::is.finite(kM) | !base::is.finite(LM))
  )
  
  res_f <- mc_summary(kF, LF, trueF$k, trueF$Linf, kB_F, LB_F, convF)
  res_m <- mc_summary(kM, LM, trueM$k, trueM$Linf, kB_M, LB_M, convM)
  
  res <- base::rbind(
    base::cbind(Sex = "Female", Scenario = label, res_f),
    base::cbind(Sex = "Male",   Scenario = label, res_m)
  )
  rownames(res) <- NULL
  
  list(summary = res, diagnostics = base::rbind(diag_f, diag_m))
}

base::message("Running Monte Carlo (R = ", R_mc, ") ...")

out_fixed     <- run_mc(R_mc, sim_fun_fixed,
                        t_obs = t_obs_fixed, sigma_eps = sigma_eps,
                        sdlog_L = sdlog_L, shape_L = shape_L, label = "Fixed")

out_lognormal <- run_mc(R_mc, sim_fun_logn,
                        t_obs = t_obs_long, sigma_eps = sigma_eps,
                        sdlog_L = sdlog_L, shape_L = shape_L, label = "Lognormal")

out_gamma     <- run_mc(R_mc, sim_fun_gamma,
                        t_obs = t_obs_long, sigma_eps = sigma_eps,
                        sdlog_L = sdlog_L, shape_L = shape_L, label = "Gamma")

cat("\n==================== MONTE-CARLO SUMMARIES ====================\n")
cat("\n--- FIXED ---\n");     print(out_fixed$summary)
cat("\n--- LOGNORMAL ---\n"); print(out_lognormal$summary)
cat("\n--- GAMMA ---\n");     print(out_gamma$summary)
