test_table5 <- function(R = 200, n_f = 50, n_m = 50) {
  library(dplyr)
  set.seed(123)
  
  # --- True parameter values (per sex) ---
  k_true_f    <- 0.275; Linf_true_f <- 186.23
  k_true_m    <- 0.247; Linf_true_m <- 228.04
  zeta_true_f <- 75.055; zeta_true_m <- 59.947
  gam_f <- 12.535; b1_f <- -2.087; b2_f <- 0.011
  gam_m <- 11.33;  b1_m <- -1.984; b2_m <- 0.010
  
  meanlog_f_Linf <- log(Linf_true_f) - 0.5 * 0.05^2
  meanlog_m_Linf <- log(Linf_true_m) - 0.5 * 0.05^2
  shape_f_Linf <- 200; scale_f_Linf <- Linf_true_f / shape_f_Linf
  shape_m_Linf <- 200; scale_m_Linf <- Linf_true_m / shape_m_Linf
  
  # --- Simulator ---
  simulate_sex <- function(n, sex_code, k_true, Linf_true, zeta_true, gam, b1, b2,
                           re_type = c("fixed","lognormal","gamma"),
                           meanlog_Linf=NULL, sdlog_Linf=NULL,
                           shape_Linf=NULL, scale_Linf=NULL) {
    re_type <- match.arg(re_type)
    if (n <= 0) return(NULL)
    ni_vals <- if (sex_code == 1) 4:15 else 2:15
    ni <- sample(ni_vals, n, replace=TRUE)
    T0 <- sample(10:150, n, replace=TRUE) / 365
    out <- vector("list", sum(ni)); idx <- 1L
    for (i in seq_len(n)) {
      Linf_i <- switch(re_type,
                       "fixed"     = Linf_true,
                       "lognormal" = rlnorm(1, meanlog=meanlog_Linf, sdlog=sdlog_Linf),
                       "gamma"     = rgamma(1, shape=shape_Linf, scale=scale_Linf))
      PL <- Linf_i * (1 - exp(-k_true * T0[i]))
      Tcum <- T0[i]
      for (j in seq_len(ni[i])) {
        a <- (1 - exp(-k_true * Tcum)) * (zeta_true - 1)
        b <- exp(-k_true * Tcum) * (zeta_true - 1)
        if (a <= 0 || b <= 0) next
        Lambda <- rbeta(1, a, b)
        INC <- Lambda * (Linf_i - PL)
        mu_IP <- exp(b1 + b2 * PL)
        if (mu_IP <= 0) next
        IP <- rgamma(1, shape=gam*mu_IP, scale=1/gam)
        out[[idx]] <- data.frame(LOBSTER=i, SEX=sex_code, PL=PL, INC=INC,
                                 INT=IP*365.25, Tcum=Tcum)
        idx <- idx+1L
        PL <- PL+INC; Tcum <- Tcum+IP
      }
    }
    bind_rows(out)
  }
  
  # --- Log-likelihood in Linf ---
  loglik_Linf <- function(Linf, dat, k_true, zeta_true, gam, b1, b2) {
    if (Linf <= max(dat$PL)) return(-Inf)
    LL <- 0
    for (iid in unique(dat$LOBSTER)) {
      dati <- dat[dat$LOBSTER==iid,]
      Xi <- dati$INC; Li <- dati$PL; IP_i <- dati$INT/365.25; Tcum <- dati$Tcum
      for (j in seq_along(Xi)) {
        a <- (1 - exp(-k_true*Tcum[j]))*(zeta_true-1)
        b <- exp(-k_true*Tcum[j])*(zeta_true-1)
        if (a <= 0 || b <= 0) return(-Inf)
        denom <- Linf - Li[j]; if (denom <= 0) return(-Inf)
        x <- Xi[j]/denom; if (x <= 0 || x >= 1) return(-Inf)
        log_beta <- dbeta(x,a,b,log=TRUE) - log(denom)
        mu_IP <- exp(b1 + b2*Li[j]); shape <- gam*mu_IP; scale <- 1/gam
        if (shape <= 0) return(-Inf)
        log_gamma <- dgamma(IP_i[j], shape=shape, scale=scale, log=TRUE)
        LL <- LL + log_beta + log_gamma
      }
    }
    LL
  }
  
  estimate_Linf <- function(dat, k_true, zeta_true, gam, b1, b2, true_Linf) {
    if (is.null(dat) || nrow(dat)==0) return(NA_real_)
    lower <- max(true_Linf-30, max(dat$PL)+1e-3)
    upper <- true_Linf+30
    opt <- try(optimize(loglik_Linf, interval=c(lower,upper),
                        dat=dat, k_true=k_true, zeta_true=zeta_true,
                        gam=gam, b1=b1, b2=b2, maximum=TRUE), silent=TRUE)
    if (inherits(opt,"try-error") || !is.finite(opt$objective)) return(NA_real_)
    opt$maximum
  }
  
  run_once_design <- function(design, n_f=50, n_m=50) {
    dat_f <- simulate_sex(n_f,1,k_true_f,Linf_true_f,zeta_true_f,gam_f,b1_f,b2_f,
                          re_type=design, meanlog_Linf=meanlog_f_Linf, sdlog_Linf=0.05,
                          shape_Linf=shape_f_Linf, scale_Linf=scale_f_Linf)
    Linf_hat_f <- estimate_Linf(dat_f,k_true_f,zeta_true_f,gam_f,b1_f,b2_f,Linf_true_f)
    dat_m <- simulate_sex(n_m,2,k_true_m,Linf_true_m,zeta_true_m,gam_m,b1_m,b2_m,
                          re_type=design, meanlog_Linf=meanlog_m_Linf, sdlog_Linf=0.05,
                          shape_Linf=shape_m_Linf, scale_Linf=scale_m_Linf)
    Linf_hat_m <- estimate_Linf(dat_m,k_true_m,zeta_true_m,gam_m,b1_m,b2_m,Linf_true_m)
    c(Linf_f=Linf_hat_f,Linf_m=Linf_hat_m)
  }
  
  summarise_param <- function(est,true) {
    est <- est[is.finite(est)]
    if (length(est)<3) return(c(mean=NA,bias=NA,rmse=NA))
    c(mean=mean(est), bias=mean(est-true), rmse=sqrt(mean((est-true)^2)))
  }
  
  designs <- c("fixed","lognormal","gamma")
  results <- data.frame(Design=character(),Sex=character(),
                        True=numeric(),Mean=numeric(),Bias=numeric(),RMSE=numeric(),
                        stringsAsFactors=FALSE)
  
  for (d in designs) {
    store_f <- rep(NA_real_,R); store_m <- rep(NA_real_,R)
    for (r in 1:R) {
      res <- run_once_design(d,n_f,n_m)
      store_f[r] <- res["Linf_f"]; store_m[r] <- res["Linf_m"]
    }
    res_f <- summarise_param(store_f,Linf_true_f)
    res_m <- summarise_param(store_m,Linf_true_m)
    results <- rbind(results,
                     data.frame(Design=d,Sex="Female",True=Linf_true_f,
                                Mean=res_f["mean"],Bias=res_f["bias"],RMSE=res_f["rmse"]),
                     data.frame(Design=d,Sex="Male",True=Linf_true_m,
                                Mean=res_m["mean"],Bias=res_m["bias"],RMSE=res_m["rmse"]))
  }
  
  if (!dir.exists("results/tables")) dir.create("results/tables",recursive=TRUE)
  write.csv(results,"results/tables/table5.csv",row.names=FALSE)
  
  invisible(results)
}
