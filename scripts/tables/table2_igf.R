test_table2_igf <- function() {
    

    dat <- lobster
    
    log_dinvgauss <- function(x, mu, lambda) {
        0.5 * (log(lambda) - log(2 * pi) - 3 * log(x)) -
            lambda * (x - mu)^2 / (2 * mu^2 * x)
    }
    
    LL_IG_fixed <- function(par, Linf, L, I, T) {
        
        lambda <- par[1]
        k      <- par[2]
        
        if (lambda <= 0 || k <= 0)
            return(1e10)
        
        mu <- (Linf - L) * (1 - exp(-k * T))
        
        if (any(mu <= 0) || any(I <= 0))
            return(1e10)
        
        ll <- log_dinvgauss(I, mu, lambda)
        
        return(-sum(ll))
    }
    
    fit_IGF_fixed <- function(df, Linf, start, lower, upper) {
        
        T <- df$NINT / 365.25
        I <- df$INC
        L <- df$CL
        
        res <- optim(
            par    = start,           # (lambda, k)
            fn     = LL_IG_fixed,
            Linf   = Linf,
            L      = L,
            I      = I,
            T      = T,
            method = "L-BFGS-B",
            lower  = lower,
            upper  = upper,
            control = list(maxit = 2000)
        )
        
        logLik <- -res$value
        p      <- 2                # FIXED Linf
        AIC    <- 2 * p - 2 * logLik
        
        list(
            par    = res$par,
            logLik = logLik,
            AIC    = AIC
        )
    }
    
    df_f <- dat[dat$SEX == 1, ]
    Linf_f <- max(df_f$CL) + 10   # explicit fixed choice
    
    fit_f <- fit_IGF_fixed(
        df    = df_f,
        Linf  = Linf_f,
        start = c(1.0, 0.05),
        lower = c(0.01, 0.01),
        upper = c(10.0, 0.5)
    )
    
    df_m <- dat[dat$SEX == 2, ]
    Linf_m <- max(df_m$CL) + 50   # explicit fixed choice
    
    fit_m <- fit_IGF_fixed(
        df    = df_m,
        Linf  = Linf_m,
        start = c(1.0, 0.05),
        lower = c(0.01, 0.01),
        upper = c(10.0, 1.0)
    )
    
    table2_results <- data.frame(
        Sex    = c("Female", "Male"),
        Linf   = c(Linf_f, Linf_m),      # FIXED values
        lambda = c(fit_f$par[1], fit_m$par[1]),
        k      = c(fit_f$par[2], fit_m$par[2]),
        logLik = c(fit_f$logLik, fit_m$logLik),
        AIC    = c(fit_f$AIC, fit_m$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        table2_results,
        "results/tables/table2_igf_fixedLinf.csv",
        row.names = FALSE
    )
    
    return(table2_results)
}

test_table2_igf()
