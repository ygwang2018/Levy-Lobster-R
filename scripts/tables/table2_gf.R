test_table2_gf <- function() {
    
    
    dat <- lobster
    
    LL_gammaMI_fixed <- function(par, Linf, L, I, T) {
        
        lambda <- par[1]   # shape-rate parameter
        k      <- par[2]
        
        if (lambda <= 0 || k <= 0)
            return(1e10)
        
        mu <- (Linf - L) * (1 - exp(-k * T))
        
        if (any(mu <= 0) || any(I <= 0))
            return(1e10)
        
        shape <- lambda * mu
        if (any(shape <= 0))
            return(1e10)
        
        # Gamma density with rate = lambda
        ll <- shape * log(lambda) +
            (shape - 1) * log(I) -
            lambda * I -
            lgamma(shape)
        
        return(-sum(ll))
    }
    
    fit_gammaMI_fixed <- function(df, Linf, start, lower, upper) {
        
        T <- df$NINT / 365.25
        I <- df$INC
        L <- df$CL
        
        res <- optim(
            par    = start,            # (lambda, k)
            fn     = LL_gammaMI_fixed,
            Linf   = Linf,
            L      = L,
            I      = I,
            T      = T,
            method = "L-BFGS-B",
            lower  = lower,
            upper  = upper
        )
        
        logLik <- -res$value
        p      <- 2                  # FIXED Linf
        AIC    <- 2 * p - 2 * logLik
        
        list(
            par    = res$par,
            logLik = logLik,
            AIC    = AIC
        )
    }
    
    df_f <- dat[dat$SEX == 1, ]
    Linf_f <- max(df_f$CL) + 10   # fixed, explicit
    
    fit_f <- fit_gammaMI_fixed(
        df    = df_f,
        Linf  = Linf_f,
        start = c(0.5, 0.05),
        lower = c(0.1, 0.01),
        upper = c(1.5, 0.5)
    )
    
    df_m <- dat[dat$SEX == 2, ]
    Linf_m <- max(df_m$CL) + 50   # fixed, explicit
    
    fit_m <- fit_gammaMI_fixed(
        df    = df_m,
        Linf  = Linf_m,
        start = c(0.5, 0.05),
        lower = c(0.1, 0.01),
        upper = c(1.5, 0.5)
    )
    
    out <- data.frame(
        Sex    = c("Female", "Male"),
        Linf   = c(Linf_f, Linf_m),      # FIXED
        lambda = c(fit_f$par[1], fit_m$par[1]),
        k      = c(fit_f$par[2], fit_m$par[2]),
        logLik = c(fit_f$logLik, fit_m$logLik),
        AIC    = c(fit_f$AIC, fit_m$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        out,
        "results/tables/table2_gammaMI_fixedLinf.csv",
        row.names = FALSE
    )
    
    return(out)
}
test_table2_gf()
