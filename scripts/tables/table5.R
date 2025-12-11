test_table5 <- function(R=200, n_f=50, n_m=50) {
  library(dplyr)
  set.seed(123)

  kf <- 0.275; Lf <- 186.23
  km <- 0.247; Lm <- 228.04
  zf <- 75.055; zm <- 59.947
  gamf <- 12.535; b1f <- -2.087; b2f <- 0.011
  gamm <- 11.33;  b1m <- -1.984; b2m <- 0.010

  meanlogLf <- log(Lf) - 0.5*0.05^2
  meanlogLm <- log(Lm) - 0.5*0.05^2
  shLf <- 200; scLf <- Lf/shLf
  shLm <- 200; scLm <- Lm/shLm

  simulate <- function(n, sex, k, Ltrue, z, gam, b1, b2,
                       type, meanlog, sdlog, shape, scale) {

    if (n <= 0) return(NULL)
    ni <- if (sex == 1) sample(3:7, n, TRUE) else sample(2:6, n, TRUE)
    T0 <- sample(10:150, n, TRUE)/365
    out <- vector("list", sum(ni)); id <- 1L

    for (i in seq_len(n)) {
      Li <- switch(type,
                   fixed=Ltrue,
                   lognormal=rlnorm(1, meanlog, sdlog),
                   gamma=rgamma(1, shape, scale))

      PL <- Li*(1 - exp(-k*T0[i]))
      Tc <- T0[i]

      for (j in seq_len(ni[i])) {
        a <- (1-exp(-k*Tc))*(z-1)
        b <-  exp(-k*Tc)*(z-1)
        if (a<=0 || b<=0) next

        lam <- rbeta(1, a, b)
        INCt <- lam*(Li - PL)
        INC <- max(INCt + rnorm(1,0,1), 0.1)

        mu <- exp(b1 + b2*PL)
        IPt <- rgamma(1, gam*mu, scale=1/gam)
        IP  <- IPt * rgamma(1, 50, scale=1/50)

        PLobs <- PL + rnorm(1,0,1.5)

        out[[id]] <- data.frame(LOBSTER=i, SEX=sex,
                                PL=PLobs, INC=INC,
                                INT=IP*365.25, Tcum=Tc)
        id <- id+1L
        PL <- PL + INCt
        Tc <- Tc + IPt
      }
    }
    bind_rows(out)
  }

  loglik <- function(Linf, dat, k, z, gam, b1, b2) {
    if (Linf <= max(dat$PL)) return(-Inf)
    LL <- 0
    ids <- unique(dat$LOBSTER)

    for (i in ids) {
      d <- dat[dat$LOBSTER==i,]
      Xi <- d$INC; Li <- d$PL
      IPi <- d$INT/365.25
      Tc <- d$Tcum

      for (j in seq_along(Xi)) {
        a <- (1-exp(-k*Tc[j]))*(z-1)
        b <-  exp(-k*Tc[j])*(z-1)
        if (a<=0 || b<=0) return(-Inf)

        den <- Linf - Li[j]
        if (den<=0) return(-Inf)
        x <- Xi[j]/den
        if (x<=0 || x>=1) return(-Inf)

        logb <- dbeta(x, a, b, log=TRUE) - log(den)
        mu <- exp(b1 + b2*Li[j])
        sh <- gam*mu
        lg <- dgamma(IPi[j], sh, scale=1/gam, log=TRUE)
        LL <- LL + logb + lg
      }
    }
    LL
  }

  estLinf <- function(dat, k, z, gam, b1, b2) {
    if (is.null(dat) || nrow(dat)==0) return(NA_real_)
    lo <- max(dat$PL) + 1e-3
    up <- lo + 80
    out <- try(optimize(loglik, c(lo,up), dat=dat,
                        k=k, z=z, gam=gam, b1=b1, b2=b2,
                        maximum=TRUE), silent=TRUE)
    if (inherits(out,"try-error") || !is.finite(out$objective))
      return(NA_real_)
    out$maximum
  }

  sumstat <- function(est, true) {
    est <- est[is.finite(est)]
    if (length(est)<3) return(c(mean=NA, bias=NA, var=NA, se=NA, rmse=NA))
    m  <- mean(est)
    b  <- mean(est - true)
    v  <- var(est)
    s  <- sqrt(v)
    r  <- sqrt(mean((est - true)^2))
    c(mean=m, bias=b, var=v, se=s, rmse=r)
  }

  designs <- c("fixed","lognormal","gamma")
  out <- data.frame()

  for (d in designs) {
    Fhat <- numeric(R)
    Mhat <- numeric(R)

    for (r in 1:R) {
      df <- simulate(n_f,1,kf,Lf,zf,gamf,b1f,b2f,d,
                     meanlogLf,0.05,shLf,scLf)
      dm <- simulate(n_m,2,km,Lm,zm,gamm,b1m,b2m,d,
                     meanlogLm,0.05,shLm,scLm)

      Fhat[r] <- estLinf(df,kf,zf,gamf,b1f,b2f)
      Mhat[r] <- estLinf(dm,km,zm,gamm,b1m,b2m)
    }

    sf <- sumstat(Fhat, Lf)
    sm <- sumstat(Mhat, Lm)

    out <- rbind(out,
      data.frame(Design=d, Sex="Female", True=Lf,
                 Mean=sf["mean"], Bias=sf["bias"],
                 Variance=sf["var"], SE=sf["se"], RMSE=sf["rmse"]),
      data.frame(Design=d, Sex="Male", True=Lm,
                 Mean=sm["mean"], Bias=sm["bias"],
                 Variance=sm["var"], SE=sm["se"], RMSE=sm["rmse"])
    )
  }

  if (!dir.exists("results/tables")) dir.create("results/tables", TRUE)
  write.csv(out, "results/tables/table5.csv", row.names=FALSE)

  invisible(out)
}
