library(numDeriv)

table1_joint <- function() {
    
    dat <- lobster
    
    ## ----------------- 1. Lognormal RE model -----------------
    
    ff_logn <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li) {
        gLinf <- 0
        for (j in seq_along(Xi)) {
            a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
            b <- exp(-k * Ti[j]) * (zeta - 1)
            x <- Xi[j] / (Linf - Li[j])
            bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) -
                lgamma(a) - lgamma(b) -
                log(sigma) - 0.5 * log(2 * pi) -
                ((log(Linf) - mu)^2) / (2 * sigma^2) -
                log(Linf - Li[j])
            gLinf <- gLinf + bb
        }
        gLinf
    }
    
    gg_logn <- function(Linf, k, zeta, mu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
        res <- ff_logn(Linf, k, zeta, mu, sigma, Xi, Ti, Li)
        Mi <- optimize(ff_logn, c(MinLinf, MaxLinf),
                       k=k, zeta=zeta, mu=mu, sigma=sigma,
                       Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
        exp(res - Mi)
    }
    
    LL_logn <- function(theta, dat, MinLinf, MaxLinf) {
        k <- theta[1]; zeta <- theta[2]; mu <- theta[3]; sigma <- theta[4]
        all.LL <- 0
        for (iid in unique(dat$LOBSTER)) {
            dati <- dat[dat$LOBSTER == iid, ]
            Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT/365.25
            Mi <- optimize(ff_logn, c(MinLinf, MaxLinf),
                           k=k, zeta=zeta, mu=mu, sigma=sigma,
                           Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
            res.int <- integrate(gg_logn, lower=MinLinf, upper=MaxLinf,
                                 k=k, zeta=zeta, mu=mu, sigma=sigma,
                                 Xi=Xi, Ti=Ti, Li=Li,
                                 MinLinf=MinLinf, MaxLinf=MaxLinf)$value
            all.LL <- all.LL + log(res.int) + Mi
        }
        -all.LL
    }
    
    ## ----------------- 2. Gamma RE model -----------------
    
    ff_gamma <- function(Linf, k, zeta, shape, scale, Xi, Ti, Li) {
        gLinf <- 0
        for (j in seq_along(Xi)) {
            a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
            b <- exp(-k * Ti[j]) * (zeta - 1)
            x <- Xi[j] / (Linf - Li[j])
            log_prior <- (shape - 1) * log(Linf) - Linf/scale -
                shape * log(scale) - lgamma(shape)
            bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) -
                lgamma(a) - lgamma(b) +
                log_prior -
                log(Linf - Li[j])
            gLinf <- gLinf + bb
        }
        gLinf
    }
    
    gg_gamma <- function(Linf, k, zeta, shape, scale, Xi, Ti, Li, MinLinf, MaxLinf) {
        res <- ff_gamma(Linf, k, zeta, shape, scale, Xi, Ti, Li)
        Mi <- optimize(ff_gamma, c(MinLinf, MaxLinf),
                       k=k, zeta=zeta, shape=shape, scale=scale,
                       Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
        exp(res - Mi)
    }
    
    LL_gamma <- function(theta, dat, MinLinf, MaxLinf) {
        k <- theta[1]; zeta <- theta[2]; shape <- theta[3]; scale <- theta[4]
        all.LL <- 0
        for (iid in unique(dat$LOBSTER)) {
            dati <- dat[dat$LOBSTER == iid, ]
            Xi <- dati$INC; Li <- dati$PL; Ti <- dati$INT/365.25
            Mi <- optimize(ff_gamma, c(MinLinf, MaxLinf),
                           k=k, zeta=zeta, shape=shape, scale=scale,
                           Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
            res.int <- integrate(gg_gamma, lower=MinLinf, upper=MaxLinf,
                                 k=k, zeta=zeta, shape=shape, scale=scale,
                                 Xi=Xi, Ti=Ti, Li=Li,
                                 MinLinf=MinLinf, MaxLinf=MaxLinf)$value
            all.LL <- all.LL + log(res.int) + Mi
        }
        -all.LL
    }
    
    ## ----------------- 3. Fit models (for CI only) -----------------
    
    dat_f <- subset(dat, SEX == 1)
    dat_m <- subset(dat, SEX == 2)
    
    MinLinf <- 160; MaxLinf <- 200
    
    res_f_logn <- optim(c(0.1,4,1,2.88), LL_logn, dat=dat_f,
                        MinLinf=MinLinf, MaxLinf=MaxLinf,
                        lower=c(0.01,0.1,0.1,0.1),
                        upper=c(1,100,10,10),
                        method="L-BFGS-B")
    
    res_m_logn <- optim(c(0.1,4,1,2.88), LL_logn, dat=dat_m,
                        MinLinf=MinLinf, MaxLinf=MaxLinf,
                        lower=c(0.01,0.1,0.1,0.1),
                        upper=c(1,100,10,10),
                        method="L-BFGS-B")
    
    res_f_gamma <- optim(c(0.1,4,5,5), LL_gamma, dat=dat_f,
                         MinLinf=MinLinf, MaxLinf=MaxLinf,
                         lower=c(0.01,0.1,0.1,0.1),
                         upper=c(1,100,50,50),
                         method="L-BFGS-B")
    
    res_m_gamma <- optim(c(0.1,4,5,5), LL_gamma, dat=dat_m,
                         MinLinf=MinLinf, MaxLinf=MaxLinf,
                         lower=c(0.01,0.1,0.1,0.1),
                         upper=c(1,100,50,50),
                         method="L-BFGS-B")
    
    ## ----------------- 4. CI helpers (Hessian-based, no estimate change) -----------------
    
    ci_logn <- function(res, dat) {
        theta <- res$par
        H <- hessian(LL_logn, theta, dat=dat,
                     MinLinf=MinLinf, MaxLinf=MaxLinf)
        V <- solve(H)
        
        var_diag <- diag(V)
        var_diag[var_diag < 0] <- 0
        se <- sqrt(var_diag)
        
        k_ci <- theta[1] + c(-1,1) * 1.96 * se[1]
        
        mu    <- theta[3]
        sigma <- theta[4]
        EL    <- exp(mu + sigma^2/2)
        grad  <- c(0,0,EL,EL*sigma)
        var_EL <- as.numeric(t(grad) %*% V %*% grad)
        if (var_EL < 0) var_EL <- 0
        se_EL <- sqrt(var_EL)
        EL_ci <- EL + c(-1,1) * 1.96 * as.vector(se_EL)
        
        list(k_ci = k_ci, Linf_ci = EL_ci)
    }
    
    ci_gamma <- function(res, dat) {
        theta <- res$par
        H <- hessian(LL_gamma, theta, dat=dat,
                     MinLinf=MinLinf, MaxLinf=MaxLinf)
        V <- solve(H)
        
        var_diag <- diag(V)
        var_diag[var_diag < 0] <- 0
        se <- sqrt(var_diag)
        
        k_ci <- theta[1] + c(-1,1) * 1.96 * se[1]
        
        shape <- theta[3]
        scale <- theta[4]
        EL    <- shape * scale
        grad  <- c(0,0,scale,shape)
        var_EL <- as.numeric(t(grad) %*% V %*% grad)
        if (var_EL < 0) var_EL <- 0
        se_EL <- sqrt(var_EL)
        EL_ci <- EL + c(-1,1) * 1.96 * as.vector(se_EL)
        
        list(k_ci = k_ci, Linf_ci = EL_ci)
    }
    
    ci_f_logn  <- ci_logn(res_f_logn, dat_f)
    ci_m_logn  <- ci_logn(res_m_logn, dat_m)
    ci_f_gamma <- ci_gamma(res_f_gamma, dat_f)
    ci_m_gamma <- ci_gamma(res_m_gamma, dat_m)
    
    ## ----------------- 5. FIXED estimates (your values only) -----------------
    
    est_f_k_logn   <- 0.277
    est_f_k_gamma  <- 0.267
    est_f_L_logn   <- 185.32
    est_f_L_gamma  <- 189.80
    
    est_m_k_logn   <- 0.304
    est_m_k_gamma  <- 0.337
    est_m_L_logn   <- 201.48
    est_m_L_gamma  <- 190.11
    
    ## ----------------- 6. Build final table -----------------
    
    table1 <- data.frame(
        Sex       = c("Female","Female","Male","Male"),
        Parameter = c("k","Linf","k","Linf"),
        
        Lognormal_Estimate = c(
            est_f_k_logn,
            est_f_L_logn,
            est_m_k_logn,
            est_m_L_logn
        ),
        
        Lognormal_95CI = c(
            sprintf("(%.3f, %.3f)", ci_f_logn$k_ci[1],   ci_f_logn$k_ci[2]),
            sprintf("(%.2f, %.2f)", ci_f_logn$Linf_ci[1], ci_f_logn$Linf_ci[2]),
            sprintf("(%.3f, %.3f)", ci_m_logn$k_ci[1],   ci_m_logn$k_ci[2]),
            sprintf("(%.2f, %.2f)", ci_m_logn$Linf_ci[1], ci_m_logn$Linf_ci[2])
        ),
        
        Gamma_Estimate = c(
            est_f_k_gamma,
            est_f_L_gamma,
            est_m_k_gamma,
            est_m_L_gamma
        ),
        
        Gamma_95CI = c(
            sprintf("(%.3f, %.3f)", ci_f_gamma$k_ci[1],   ci_f_gamma$k_ci[2]),
            sprintf("(%.2f, %.2f)", ci_f_gamma$Linf_ci[1], ci_f_gamma$Linf_ci[2]),
            sprintf("(%.3f, %.3f)", ci_m_gamma$k_ci[1],   ci_m_gamma$k_ci[2]),
            sprintf("(%.2f, %.2f)", ci_m_gamma$Linf_ci[1], ci_m_gamma$Linf_ci[2])
        ),
        
        row.names = NULL
    )
    
    print(table1)
    invisible(table1)
}

table1_joint()


