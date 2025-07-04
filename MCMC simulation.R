set.seed(123) # For reproducibility

n <- 100
years <- 10
par(mfrow=c(2,2))

# Helper function for simulation
simulate_growth <- function(n, ni_range, T_range, Linf, k, zeta, gam, b1, b2, random_Linf=FALSE) {
  ni <- sample(ni_range, n, replace=TRUE)
  T <- sample(T_range, n, replace=TRUE) / 365
  if (random_Linf) {
    Linf_vec <- rgamma(n, shape=200, scale=1/1.155)
  } else {
    Linf_vec <- rep(Linf, n)
  }
  L <- Linf_vec * (1 - exp(-k * T))
  size <- list()
  time <- list()
  for (i in 1:n) {
    PL <- numeric(ni[i]+1)
    INT <- numeric(ni[i]+1)
    PL[1] <- L[i]
    INT[1] <- T[i]
    for (j in 1:ni[i]) {
      shape1 <- (1 - exp(-k * INT[j])) * (zeta - 1)
      shape2 <- exp(-k * INT[j]) * (zeta - 1)
      PL[j+1] <- PL[j] + (Linf_vec[i] - PL[j]) * rbeta(1, shape1, shape2)
      INT[j+1] <- INT[j] + rgamma(1, shape=gam * exp(b1 + b2 * PL[j]), scale=1/gam)
    }
    size[[i]] <- PL[1:ni[i]]
    time[[i]] <- INT[1:ni[i]]
  }
  return(list(time=time, size=size))
}

# Panel 1: Fixed Linf, Females
Linf_f <- 183.87; k_f <- 0.283; zeta_f <- 75.055; gam_f <- 12.535; b1_f <- -2.087; b2_f <- 0.011
res_f_fixed <- simulate_growth(n, 4:15, 10:150, Linf_f, k_f, zeta_f, gam_f, b1_f, b2_f, random_Linf=FALSE)
plot(NA, xlim=c(0,years), ylim=c(0,230), xlab="Years", ylab="Lengths(mm)", main="Beta process for females")
for (i in 1:n) lines(res_f_fixed$time[[i]], res_f_fixed$size[[i]], type="s", col="black")
# Mean curve
all_times <- unlist(res_f_fixed$time)
all_sizes <- unlist(res_f_fixed$size)
spline_fit <- smooth.spline(all_times, all_sizes, spar=1)
lines(spline_fit, col="cyan", lwd=3)

# Panel 2: Fixed Linf, Males
Linf_m <- 228.04; k_m <- 0.2473; zeta_m <- 59.947; gam_m <- 11.33; b1_m <- -1.984; b2_m <- 0.01
res_m_fixed <- simulate_growth(n, 2:15, 5:150, Linf_m, k_m, zeta_m, gam_m, b1_m, b2_m, random_Linf=FALSE)
plot(NA, xlim=c(0,years), ylim=c(0,230), xlab="Years", ylab="Lengths(mm)", main="Beta process for males")
for (i in 1:n) lines(res_m_fixed$time[[i]], res_m_fixed$size[[i]], type="s", col="black")
all_times <- unlist(res_m_fixed$time)
all_sizes <- unlist(res_m_fixed$size)
spline_fit <- smooth.spline(all_times, all_sizes, spar=1)
lines(spline_fit, col="cyan", lwd=3)

# Panel 3: Random Linf, Females
res_f_random <- simulate_growth(n, 4:15, 10:150, Linf_f, k_f, zeta_f, gam_f, b1_f, b2_f, random_Linf=TRUE)
plot(NA, xlim=c(0,years), ylim=c(0,230), xlab="Years", ylab="Lengths(mm)", main="Beta process for females")
for (i in 1:n) lines(res_f_random$time[[i]], res_f_random$size[[i]], type="s", col="black")
all_times <- unlist(res_f_random$time)
all_sizes <- unlist(res_f_random$size)
spline_fit <- smooth.spline(all_times, all_sizes, spar=1)
lines(spline_fit, col="cyan", lwd=3)

# Panel 4: Random Linf, Males
res_m_random <- simulate_growth(n, 2:15, 5:150, Linf_m, k_m, zeta_m, gam_m, b1_m, b2_m, random_Linf=TRUE)
plot(NA, xlim=c(0,years), ylim=c(0,230), xlab="Years", ylab="Lengths(mm)", main="Beta process for males")
for (i in 1:n) lines(res_m_random$time[[i]], res_m_random$size[[i]], type="s", col="black")
all_times <- unlist(res_m_random$time)
all_sizes <- unlist(res_m_random$size)
spline_fit <- smooth.spline(all_times, all_sizes, spar=1)
lines(spline_fit, col="cyan", lwd=3)
