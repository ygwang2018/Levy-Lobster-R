library(dplyr)

test_table2_BF <- function() {
    
    dat <- lobster
    
    negLL_BF <- function(par, T, I, L) {
        
        Linf <- par[1]
        k    <- par[2]
        zeta <- par[3]
        
        if (k <= 0 || zeta <= 1 || Linf <= max(L) + 1e-6)
            return(1e10)
        
        ll <- 0
        
        for (i in seq_along(I)) {
            
            if (I[i] >= Linf - L[i]) return(1e10)
            
            alpha <- (1 - exp(-k * T[i])) * (zeta - 1)
            beta  <- exp(-k * T[i]) * (zeta - 1)
            
            if (alpha <= 0 || beta <= 0) return(1e10)
            
            x <- I[i] / (Linf - L[i])
            if (x <= 0 || x >= 1) return(1e10)
            
            ll <- ll +
                lgamma(alpha + beta) -
                lgamma(alpha) -
                lgamma(beta) +
                (alpha - 1) * log(x) +
                (beta  - 1) * log(1 - x) -
                log(Linf - L[i])
        }
        
        return(-ll)
    }
    
    fit_BF <- function(df) {
        
        T <- df$INT / 365.25
        I <- df$INC
        L <- df$PL
        
        start <- c(max(L) + 20, 0.15, 20)
        
        res <- optim(
            par    = start,
            fn     = negLL_BF,
            T      = T,
            I      = I,
            L      = L,
            method = "L-BFGS-B",
            lower  = c(max(L) + 1, 1e-4, 1 + 1e-4),
            upper  = c(max(L) + 200, 2, 200)
        )
        
        logLik <- -res$value
        p      <- 3                     # Linf, k, zeta
        AIC    <- 2 * p - 2 * logLik
        
        list(
            par    = res$par,
            logLik = logLik,
            AIC    = AIC
        )
    }
    
    fit_f <- fit_BF(dat %>% filter(SEX == 1))
    
    fit_m <- fit_BF(dat %>% filter(SEX == 2))
    
    table2 <- data.frame(
        Sex    = c("Female", "Male"),
        Linf   = c(fit_f$par[1], fit_m$par[1]),
        k      = c(fit_f$par[2], fit_m$par[2]),
        zeta   = c(fit_f$par[3], fit_m$par[3]),
        logLik = c(fit_f$logLik, fit_m$logLik),
        AIC    = c(fit_f$AIC, fit_m$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        table2,
        "results/tables/table2_BF_AIC.csv",
        row.names = FALSE
    )
    
    return(table2)
}

## Run
test_table2_BF()
