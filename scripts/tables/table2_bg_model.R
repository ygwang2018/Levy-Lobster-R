test_table2_bg <- function() {
    
    dat <- lobster
    
    # 1. Log-posterior for BG model
    ff_BG <- function(Linf, k, zeta, alpha, beta, Xi, Ti, Li) {
        
        if (Linf <= max(Li) + 1e-8) return(-1e12)
        
        # Gamma prior on Linf
        log_prior <- (alpha - 1) * log(Linf) - (Linf / beta) -
            alpha * log(beta) - lgamma(alpha)
        
        out <- log_prior
        
        for (j in seq_along(Xi)) {
            
            denom <- Linf - Li[j]
            if (denom <= 1e-12) return(-1e12)
            
            x <- Xi[j] / denom
            if (x <= 1e-12 || x >= 1 - 1e-12) return(-1e12)
            
            a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
            b <- exp(-k * Ti[j]) * (zeta - 1)
            if (a <= 0 || b <= 0) return(-1e12)
            
            out <- out + dbeta(x, a, b, log = TRUE) - log(denom)
        }
        
        return(out)
    }
    
    kernel_BG <- function(Linf, Mi, k, zeta, alpha, beta, Xi, Ti, Li) {
        exp(ff_BG(Linf, k, zeta, alpha, beta, Xi, Ti, Li) - Mi)
    }
    
    safe_integrate <- function(fun, lower, upper, n = 20) {
        xs <- seq(lower, upper, length.out = n)
        vals <- sapply(xs, fun)
        dx <- (upper - lower) / (n - 1)
        sum(vals) * dx
    }
    
    # 2. Marginal log-likelihood
    LL_BG <- function(theta, dat, MinLinf, MaxLinf) {
        
        k     <- theta[1]
        zeta  <- theta[2]
        alpha <- theta[3]
        beta  <- theta[4]
        
        if (k <= 0 || zeta <= 1.01 || alpha <= 0 || beta <= 0)
            return(1e12)
        
        LLtot <- 0
        
        for (id in unique(dat$LOBSTER)) {
            
            dati <- dat[dat$LOBSTER == id, ]
            Xi <- dati$INC
            Li <- dati$PL
            Ti <- dati$INT / 365.25
            
            Mi <- optimize(
                function(Ls) ff_BG(Ls, k, zeta, alpha, beta, Xi, Ti, Li),
                interval = c(MinLinf, MaxLinf),
                maximum = TRUE
            )$objective
            
            val <- safe_integrate(
                function(Ls) kernel_BG(Ls, Mi, k, zeta, alpha, beta, Xi, Ti, Li),
                MinLinf, MaxLinf
            )
            
            if (!is.finite(val) || val <= 0) return(1e12)
            
            LLtot <- LLtot + log(val) + Mi
        }
        
        return(-LLtot)
    }
    
    # 3. SEs for k and Linf only (delta method)
    get_se_k_Linf <- function(opt, dat, MinLinf, MaxLinf) {
        
        hess <- optimHess(
            opt$par, fn = LL_BG,
            dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf
        )
        
        vcov <- tryCatch(
            solve(hess),
            error = function(e) matrix(NA, nrow = 4, ncol = 4)
        )
        
        # SE for k
        se_k <- sqrt(vcov[1, 1])
        
        # Delta-method SE for Linf = alpha * beta
        alpha <- opt$par[3]
        beta  <- opt$par[4]
        
        var_Linf <- beta^2 * vcov[3, 3] +
            alpha^2 * vcov[4, 4] +
            2 * alpha * beta * vcov[3, 4]
        
        se_Linf <- sqrt(var_Linf)
        
        list(se_k = se_k, se_Linf = se_Linf)
    }
    
    # 4. FEMALE
    dat_female <- dat[dat$SEX == 1, ]
    MinLinf <- 180
    MaxLinf <- 200
    
    init_theta_f <- c(k = 0.1, zeta = 4, alpha = 170, beta = 1.0)
    
    res_female_BG <- optim(
        par     = init_theta_f,
        fn      = LL_BG,
        dat     = dat_female,
        MinLinf = MinLinf,
        MaxLinf = MaxLinf,
        method  = "L-BFGS-B",
        lower   = c(0.01, 1.02, 1e-4, 1e-4),
        upper   = c(2, 200, 500, 20),
        control = list(maxit = 1500)
    )
    
    val_female <- res_female_BG$par
    Linf_female <- val_female[3] * val_female[4]
    se_female <- get_se_k_Linf(res_female_BG, dat_female, MinLinf, MaxLinf)
    
    AIC_female <- 2 * res_female_BG$value + 2 * length(val_female)
    
    # 5. MALE
    dat_male <- dat[dat$SEX == 2, ]
    
    init_theta_m <- c(k = 0.1, zeta = 4, alpha = 10, beta = 17)
    
    res_male_BG <- optim(
        par     = init_theta_m,
        fn      = LL_BG,
        dat     = dat_male,
        MinLinf = MinLinf,
        MaxLinf = MaxLinf,
        method  = "L-BFGS-B",
        lower   = c(0.01, 1.02, 1, 1),
        upper   = c(1, 100, 500, 50),
        control = list(maxit = 5000)
    )
    
    val_male <- res_male_BG$par
    Linf_male <- val_male[3] * val_male[4]
    se_male <- get_se_k_Linf(res_male_BG, dat_male, MinLinf, MaxLinf)
    
    AIC_male <- 2 * res_male_BG$value + 2 * length(val_male)
    
    # 6. FINAL TABLE (ONLY k & Linf SEs)
    table2_results <- data.frame(
        Sex      = c("Female", "Male"),
        k        = c(val_female[1], val_male[1]),
        k_SE     = c(se_female$se_k, se_male$se_k),
        Linf     = c(Linf_female, Linf_male),
        Linf_SE  = c(se_female$se_Linf, se_male$se_Linf),
        AIC      = c(AIC_female, AIC_male)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive = TRUE)
    
    write.csv(
        table2_results,
        "results/tables/table2_bg_model.csv",
        row.names = FALSE
    )
    
    print(table2_results)
    invisible(table2_results)
}
test_table2_bg()
