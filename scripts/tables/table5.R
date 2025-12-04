library(dplyr)

set.seed(123)

## 1. True parameter values (per sex)

k_true_f    <- 0.275
Linf_true_f <- 186.23

k_true_m    <- 0.247
Linf_true_m <- 228.04

## Beta-subordinator zeta values from your fitted model
zeta_true_f <- 75.055
zeta_true_m <- 59.947

## Gamma–GLM IP parameters from your fitted model
gam_f  <- 12.535
b1_f   <- -2.087
b2_f   <-  0.011

gam_m  <- 11.33
b1_m   <- -1.984
b2_m   <-  0.010

## Random-effects parameters
## For lognormal L∞:
sdlog_f_Linf <- 0.05
sdlog_m_Linf <- 0.05

meanlog_f_Linf <- log(Linf_true_f) - 0.5 * sdlog_f_Linf^2
meanlog_m_Linf <- log(Linf_true_m) - 0.5 * sdlog_m_Linf^2

## For gamma L∞:
shape_f_Linf <- 200
scale_f_Linf <- Linf_true_f / shape_f_Linf

shape_m_Linf <- 200
scale_m_Linf <- Linf_true_m / shape_m_Linf


## 2. Simulator: MI–IP process with optional random-effects in L∞
simulate_sex <- function(n, sex_code,
                         k_true, Linf_true,
                         zeta_true, gam, b1, b2,
                         re_type = c("fixed", "lognormal", "gamma"),
                         meanlog_Linf = NULL, sdlog_Linf = NULL,
                         shape_Linf   = NULL, scale_Linf = NULL) {
    
    re_type <- match.arg(re_type)
    
    if (n <= 0) return(NULL)
    
    ni_vals <- if (sex_code == 1) 4:15 else 2:15
    ni      <- sample(ni_vals, n, replace = TRUE)
    
    # initial time T0 in years
    T0 <- sample(10:150, n, replace = TRUE) / 365
    
    out <- vector("list", sum(ni))
    idx <- 1L
    
    for (i in seq_len(n)) {
        
        ## Choose individual L∞i according to design
        Linf_i <-
            switch(re_type,
                   "fixed"     = Linf_true,
                   "lognormal" = rlnorm(1, meanlog = meanlog_Linf, sdlog = sdlog_Linf),
                   "gamma"     = rgamma(1, shape = shape_Linf, scale = scale_Linf))
        
        PL   <- Linf_i * (1 - exp(-k_true * T0[i]))  # starting premoult length
        Tcum <- T0[i]                                 # cumulative time in years
        
        for (j in seq_len(ni[i])) {
            
            ## 1) Beta-subordinator increment
            a <- (1 - exp(-k_true * Tcum)) * (zeta_true - 1)
            b <- exp(-k_true * Tcum)       * (zeta_true - 1)
            if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) next
            
            Lambda <- rbeta(1, a, b)
            INC    <- Lambda * (Linf_i - PL)
            
            ## 2) Gamma–GLM IP
            mu_IP <- exp(b1 + b2 * PL)
            if (!is.finite(mu_IP) || mu_IP <= 0) next
            
            IP <- rgamma(1, shape = gam * mu_IP, scale = 1/gam)
            
            out[[idx]] <- data.frame(
                LOBSTER = i,
                SEX     = sex_code,
                PL      = PL,
                INC     = INC,
                INT     = IP * 365.25,   # store IP in days
                Tcum    = Tcum           # cumulative time in years
            )
            idx <- idx + 1L
            
            PL   <- PL + INC
            Tcum <- Tcum + IP
        }
    }
    
    bind_rows(out)
}

## 3. Log-likelihood in L∞ with k, zeta, (gam, b1, b2) fixed
loglik_Linf <- function(Linf, dat, k_true, zeta_true, gam, b1, b2) {
    
    if (!is.finite(Linf) || Linf <= max(dat$PL)) return(-Inf)
    
    LL <- 0
    
    for (iid in unique(dat$LOBSTER)) {
        dati <- dat[dat$LOBSTER == iid, ]
        
        Xi   <- dati$INC
        Li   <- dati$PL
        IP_i <- dati$INT / 365.25    # years
        Tcum <- dati$Tcum            # years
        
        for (j in seq_along(Xi)) {
            
            ## 1) Beta MI term
            a <- (1 - exp(-k_true * Tcum[j])) * (zeta_true - 1)
            b <- exp(-k_true * Tcum[j])       * (zeta_true - 1)
            if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(-Inf)
            
            denom <- Linf - Li[j]
            if (!is.finite(denom) || denom <= 0) return(-Inf)
            
            x <- Xi[j] / denom
            if (!is.finite(x) || x <= 0 || x >= 1) return(-Inf)
            
            log_beta <- dbeta(x, a, b, log = TRUE) - log(denom)
            if (!is.finite(log_beta)) return(-Inf)
            
            ## 2) Gamma IP term
            mu_IP <- exp(b1 + b2 * Li[j])
            shape <- gam * mu_IP
            scale <- 1 / gam
            if (!is.finite(shape) || shape <= 0) return(-Inf)
            
            log_gamma <- dgamma(IP_i[j], shape = shape, scale = scale, log = TRUE)
            if (!is.finite(log_gamma)) return(-Inf)
            
            LL <- LL + log_beta + log_gamma
        }
    }
    
    LL
}

## 4. Estimate population L∞ by 1D optimisation
estimate_Linf <- function(dat, k_true, zeta_true, gam, b1, b2, true_Linf) {
    
    if (is.null(dat) || nrow(dat) == 0) return(NA_real_)
    
    # base interval around true value (simulation knows the truth)
    lower <- true_Linf - 30
    upper <- true_Linf + 30
    
    # ensure L∞ > max(PL) + epsilon
    maxPL <- max(dat$PL)
    eps   <- 1e-3
    if (lower <= maxPL + eps) {
        lower <- maxPL + eps
    }
    
    opt <- try(
        optimize(
            f        = loglik_Linf,
            interval = c(lower, upper),
            dat      = dat,
            k_true   = k_true,
            zeta_true = zeta_true,
            gam      = gam,
            b1       = b1,
            b2       = b2,
            maximum  = TRUE
        ),
        silent = TRUE
    )
    
    if (inherits(opt, "try-error")) return(NA_real_)
    if (!is.finite(opt$objective))  return(NA_real_)
    
    opt$maximum
}

## 5. One Monte Carlo replication for a given design
run_once_design <- function(design, n_f = 50, n_m = 50) {
    
    design <- match.arg(design, c("fixed", "lognormal", "gamma"))
    
    ## Females
    dat_f <- switch(
        design,
        "fixed" = simulate_sex(
            n = n_f, sex_code = 1,
            k_true = k_true_f, Linf_true = Linf_true_f,
            zeta_true = zeta_true_f, gam = gam_f, b1 = b1_f, b2 = b2_f,
            re_type = "fixed"
        ),
        "lognormal" = simulate_sex(
            n = n_f, sex_code = 1,
            k_true = k_true_f, Linf_true = Linf_true_f,
            zeta_true = zeta_true_f, gam = gam_f, b1 = b1_f, b2 = b2_f,
            re_type = "lognormal",
            meanlog_Linf = meanlog_f_Linf, sdlog_Linf = sdlog_f_Linf
        ),
        "gamma" = simulate_sex(
            n = n_f, sex_code = 1,
            k_true = k_true_f, Linf_true = Linf_true_f,
            zeta_true = zeta_true_f, gam = gam_f, b1 = b1_f, b2 = b2_f,
            re_type = "gamma",
            shape_Linf = shape_f_Linf, scale_Linf = scale_f_Linf
        )
    )
    
    Linf_hat_f <- estimate_Linf(
        dat       = dat_f,
        k_true    = k_true_f,
        zeta_true = zeta_true_f,
        gam       = gam_f,
        b1        = b1_f,
        b2        = b2_f,
        true_Linf = Linf_true_f
    )
    
    ## Males
    dat_m <- switch(
        design,
        "fixed" = simulate_sex(
            n = n_m, sex_code = 2,
            k_true = k_true_m, Linf_true = Linf_true_m,
            zeta_true = zeta_true_m, gam = gam_m, b1 = b1_m, b2 = b2_m,
            re_type = "fixed"
        ),
        "lognormal" = simulate_sex(
            n = n_m, sex_code = 2,
            k_true = k_true_m, Linf_true = Linf_true_m,
            zeta_true = zeta_true_m, gam = gam_m, b1 = b1_m, b2 = b2_m,
            re_type = "lognormal",
            meanlog_Linf = meanlog_m_Linf, sdlog_Linf = sdlog_m_Linf
        ),
        "gamma" = simulate_sex(
            n = n_m, sex_code = 2,
            k_true = k_true_m, Linf_true = Linf_true_m,
            zeta_true = zeta_true_m, gam = gam_m, b1 = b1_m, b2 = b2_m,
            re_type = "gamma",
            shape_Linf = shape_m_Linf, scale_Linf = scale_m_Linf
        )
    )
    
    Linf_hat_m <- estimate_Linf(
        dat       = dat_m,
        k_true    = k_true_m,
        zeta_true = zeta_true_m,
        gam       = gam_m,
        b1        = b1_m,
        b2        = b2_m,
        true_Linf = Linf_true_m
    )
    
    c(Linf_f = Linf_hat_f, Linf_m = Linf_hat_m)
}

## 6. Monte Carlo over all three designs
R <- 200  # number of replications per design

designs <- c("fixed", "lognormal", "gamma")

store <- lapply(designs, function(x) list(f = rep(NA_real_, R),
                                          m = rep(NA_real_, R)))
names(store) <- designs

for (d in designs) {
    cat("\n=== Design:", d, "===\n")
    for (r in 1:R) {
        cat("  Replication", r, "of", R, "\n")
        res <- run_once_design(d, n_f = 50, n_m = 50)
        store[[d]]$f[r] <- res["Linf_f"]
        store[[d]]$m[r] <- res["Linf_m"]
    }
}

## 7. Summaries: mean, bias, RMSE for each design × sex
summarise_param <- function(est, true) {
    est <- est[is.finite(est)]
    cat("  Number of finite estimates:", length(est), "\n")
    if (length(est) < 3) {
        return(c(mean = NA, bias = NA, rmse = NA))
    }
    mean_est <- mean(est)
    bias     <- mean(est - true)
    rmse     <- sqrt(mean((est - true)^2))
    c(mean = mean_est, bias = bias, rmse = rmse)
}

results <- data.frame(
    Design = character(),
    Sex    = character(),
    True   = numeric(),
    Mean   = numeric(),
    Bias   = numeric(),
    RMSE   = numeric(),
    stringsAsFactors = FALSE
)

for (d in designs) {
    # Monte Carlo Summary
    
    cat("Female L∞\n")
    res_f <- summarise_param(store[[d]]$f, Linf_true_f)
    print(res_f)
    
    cat("Male L∞\n")
    res_m <- summarise_param(store[[d]]$m, Linf_true_m)
    print(res_m)
    
    results <- rbind(
        results,
        data.frame(
            Design = d, Sex = "Female",
            True = Linf_true_f,
            Mean = res_f["mean"],
            Bias = res_f["bias"],
            RMSE = res_f["rmse"]
        ),
        data.frame(
            Design = d, Sex = "Male",
            True = Linf_true_m,
            Mean = res_m["mean"],
            Bias = res_m["bias"],
            RMSE = res_m["rmse"]
        )
    )
}

print(results)

## 8. AOAS-style LaTeX rows
make_row <- function(design_label, sex, true_L, mean_L, bias_L, rmse_L) {
    sprintf(
        "%s & %s & %.2f & %.2f & %.2f & %.2f \\\\",
        sex, design_label, true_L, mean_L, bias_L, rmse_L
    )
}

cat("\nLaTeX rows for AOAS simulation table:\n\n")
for (i in seq_len(nrow(results))) {
    d_lab <- switch(results$Design[i],
                    "fixed"     = "Fixed $L_\\infty$",
                    "lognormal" = "Lognormal RE $L_\\infty$",
                    "gamma"     = "Gamma RE $L_\\infty$")
    cat(
        make_row(
            d_lab,
            results$Sex[i],
            results$True[i],
            results$Mean[i],
            results$Bias[i],
            results$RMSE[i]
        ),
        "\n"
    )
}


# Save results dataframe as CSV
write.csv(
    results,
    file = "results/tables/table5.csv",
    row.names = FALSE
)

