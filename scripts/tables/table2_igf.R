test_table2_igf <- function() {
    
    dat <- lobster
    
    log_dinvgauss <- function(x, mu, lambda) {
        0.5 * (log(lambda) - log(2*pi) - 3*log(x)) -
            lambda * (x - mu)^2 / (2 * mu^2 * x)
    }
    
    LL_IG <- function(par, L, I, T) {
        
        Linf   <- par[1]
        lambda <- par[2]
        k      <- par[3]
        
        if (Linf <= max(L) + 1e-6 || lambda <= 0 || k <= 0)
            return(1e10)
        
        mu <- (Linf - L) * (1 - exp(-k * T))
        
        if (any(mu <= 0) || any(I <= 0))
            return(1e10)
        
        ll <- log_dinvgauss(I, mu, lambda)
        
        return(-sum(ll))
    }
    
    fit_IGF <- function(df, start, lower, upper) {
        
        T <- df$INT / 365.25
        I <- df$INC
        L <- df$PL
        
        res <- optim(
            par     = start,
            fn      = LL_IG,
            L       = L,
            I       = I,
            T       = T,
            method  = "L-BFGS-B",
            lower   = lower,
            upper   = upper,
            control = list(maxit = 2000)
        )
        
        logLik <- -res$value
        p      <- 3
        AIC    <- 2 * p - 2 * logLik
        
        list(
            par    = res$par,
            logLik = logLik,
            AIC    = AIC
        )
    }
    
    ## Female
    fit_f <- fit_IGF(
        df    = dat[dat$SEX == 1, ],
        start = c(180, 1, 0.1),
        lower = c(100, 0.01, 0.01),
        upper = c(300, 10, 0.4)
    )
    
    ## Male
    fit_m <- fit_IGF(
        df    = dat[dat$SEX == 2, ],
        start = c(160, 1, 0.1),
        lower = c(160, 0.01, 0.01),
        upper = c(180, 10, 1)
    )
    
    table2_results <- data.frame(
        Sex    = c("Female", "Male"),
        Linf   = c(fit_f$par[1], fit_m$par[1]),
        lambda = c(fit_f$par[2], fit_m$par[2]),
        k      = c(fit_f$par[3], fit_m$par[3]),
        logLik = c(fit_f$logLik, fit_m$logLik),
        AIC    = c(fit_f$AIC, fit_m$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        table2_results,
        "results/tables/table2_igf.csv",
        row.names = FALSE
    )
    
    return(table2_results)
}
test_table2_igf()
