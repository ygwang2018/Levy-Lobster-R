test_table2_bf <- function() {
    
    
    dat <- lobster
    
    
    LL_BF_fixed <- function(par, Linf, L, I, T) {
        
        k    <- par[1]
        zeta <- par[2]
        
        if (k <= 0 || zeta <= 1)
            return(1e10)
        
        ll <- 0
        
        for (i in seq_along(I)) {
            
            denom <- Linf - L[i]
            if (denom <= 0 || I[i] >= denom)
                return(1e10)
            
            alpha <- (1 - exp(-k * T[i])) * (zeta - 1)
            beta  <- exp(-k * T[i]) * (zeta - 1)
            
            if (alpha <= 0 || beta <= 0)
                return(1e10)
            
            x <- I[i] / denom
            if (x <= 0 || x >= 1)
                return(1e10)
            
            ll <- ll +
                lgamma(alpha + beta) -
                lgamma(alpha) -
                lgamma(beta) +
                (alpha - 1) * log(x) +
                (beta  - 1) * log(1 - x) -
                log(denom)
        }
        
        return(-ll)
    }
    
    fit_BF_fixed <- function(df, Linf, start, lower, upper) {
        
        T <- df$NINT / 365.25
        I <- df$INC
        L <- df$CL
        
        res <- optim(
            par    = start,           # (k, zeta)
            fn     = LL_BF_fixed,
            Linf   = Linf,
            L      = L,
            I      = I,
            T      = T,
            method = "L-BFGS-B",
            lower  = lower,
            upper  = upper
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
    Linf_f <- max(df_f$CL) + 10   # fixed choice
    
    fit_f <- fit_BF_fixed(
        df    = df_f,
        Linf  = Linf_f,
        start = c(0.15, 20),
        lower = c(1e-4, 1 + 1e-4),
        upper = c(2.0, 200)
    )
    
    df_m <- dat[dat$SEX == 2, ]
    Linf_m <- max(df_m$CL) + 50   # fixed choice
    
    fit_m <- fit_BF_fixed(
        df    = df_m,
        Linf  = Linf_m,
        start = c(0.15, 20),
        lower = c(1e-4, 1 + 1e-4),
        upper = c(2.0, 200)
    )
    
    table2_results <- data.frame(
        Sex    = c("Female", "Male"),
        Linf   = c(Linf_f, Linf_m),      # FIXED
        k      = c(fit_f$par[1], fit_m$par[1]),
        zeta   = c(fit_f$par[2], fit_m$par[2]),
        logLik = c(fit_f$logLik, fit_m$logLik),
        AIC    = c(fit_f$AIC, fit_m$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        table2_results,
        "results/tables/table2_bf_fixedLinf.csv",
        row.names = FALSE
    )
    
    return(table2_results)
}
test_table2_bf()
