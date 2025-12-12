test_table2_gl <- function() {
    
    
    dat <- lobster
    
    ff_GL <- function(Linf, k, alpha, mu, sigma, Xi, Ti, Li) {
        
        if (!is.finite(Linf) || Linf <= max(Li) + 1e-10)
            return(-1e12)
        
        ll <- dlnorm(Linf, meanlog = mu, sdlog = sigma, log = TRUE)
        
        for (j in seq_along(Xi)) {
            mu_inc <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
            if (mu_inc <= 0 || Xi[j] <= 0)
                return(-1e12)
            ll <- ll + dgamma(Xi[j], shape = alpha,
                              scale = mu_inc / alpha, log = TRUE)
        }
        
        ll
    }
    
    kernel_GL <- function(Linf, Mi, k, alpha, mu, sigma, Xi, Ti, Li) {
        exp(ff_GL(Linf, k, alpha, mu, sigma, Xi, Ti, Li) - Mi)
    }
    
    # faster, coarse integration (sufficient for profiling)
    safe_integrate <- function(fun, lower, upper, n = 80) {
        xs <- seq(lower, upper, length.out = n)
        vals <- sapply(xs, fun)
        dx <- (upper - lower) / (n - 1)
        sum(vals) * dx
    }
    
    # ============================================================
    # 2. Marginal negative log-likelihood
    # ============================================================
    LL_GL <- function(theta, dat, MinLinf, MaxLinf) {
        
        k     <- theta[1]
        alpha <- theta[2]
        mu    <- theta[3]
        sigma <- theta[4]
        
        if (k <= 0 || alpha <= 0 || sigma <= 0)
            return(1e12)
        
        LLtot <- 0
        
        for (id in unique(dat$LOBSTER)) {
            
            dati <- dat[dat$LOBSTER == id, ]
            Xi <- dati$INC
            Li <- dati$PL
            Ti <- dati$INT / 365.25
            
            if (any(Li + Xi >= MaxLinf))
                next
            
            Mi <- optimize(
                function(Ls) ff_GL(Ls, k, alpha, mu, sigma, Xi, Ti, Li),
                interval = c(MinLinf, MaxLinf),
                maximum = TRUE
            )$objective
            
            val <- safe_integrate(
                function(Ls)
                    kernel_GL(Ls, Mi, k, alpha, mu, sigma, Xi, Ti, Li),
                MinLinf, MaxLinf
            )
            
            if (!is.finite(val) || val <= 0)
                return(1e12)
            
            LLtot <- LLtot + log(val) + Mi
        }
        
        -LLtot
    }
    
    # ============================================================
    # 3. Hessian check (for female)
    # ============================================================
    hessian_ok <- function(theta, dat, MinLinf, MaxLinf) {
        hess <- tryCatch(
            optimHess(theta, LL_GL,
                      dat = dat, MinLinf = MinLinf, MaxLinf = MaxLinf),
            error = function(e) NULL
        )
        if (is.null(hess)) return(FALSE)
        ev <- eigen(hess, symmetric = TRUE, only.values = TRUE)$values
        all(ev > 0)
    }
    
    # ============================================================
    # 4. FAST PROFILE SEs (male only)
    # ============================================================
    profile_se_k_fast <- function(fit, dat, MinLinf, MaxLinf) {
        
        k_hat <- fit$par[1]
        LLmin <- fit$value
        
        grid <- seq(0.75 * k_hat, 1.25 * k_hat, length.out = 11)
        prof <- numeric(length(grid))
        
        start <- fit$par[-1]
        
        for (i in seq_along(grid)) {
            k0 <- grid[i]
            fn <- function(rest) {
                theta <- c(k0, rest)
                LL_GL(theta, dat, MinLinf, MaxLinf)
            }
            opt <- optim(start, fn,
                         method="L-BFGS-B",
                         lower=c(0.1,0.1,0.01),
                         upper=c(100,10,10))
            prof[i] <- opt$value
            start <- opt$par
        }
        
        idx <- which(prof <= LLmin + 0.5)
        if (length(idx) < 2) return(NA_real_)
        
        (max(grid[idx]) - min(grid[idx])) / 2
    }
    
    profile_se_Linf_fast <- function(fit, dat, MinLinf, MaxLinf) {
        
        mu_hat    <- fit$par[3]
        sigma_hat <- fit$par[4]
        Linf_hat  <- exp(mu_hat + 0.5 * sigma_hat^2)
        LLmin     <- fit$value
        
        grid <- seq(Linf_hat - 25, Linf_hat + 25, length.out = 11)
        prof <- numeric(length(grid))
        
        for (i in seq_along(grid)) {
            L0 <- grid[i]
            fn <- function(sigma) {
                mu <- log(L0) - 0.5 * sigma^2
                theta <- c(fit$par[1], fit$par[2], mu, sigma)
                LL_GL(theta, dat, MinLinf, MaxLinf)
            }
            opt <- optim(sigma_hat, fn,
                         method="L-BFGS-B",
                         lower=0.01, upper=2)
            prof[i] <- opt$value
            sigma_hat <- opt$par
        }
        
        idx <- which(prof <= LLmin + 0.5)
        if (length(idx) < 2) return(NA_real_)
        
        (max(grid[idx]) - min(grid[idx])) / 2
    }
    
    # ============================================================
    # 5. Fit by sex
    # ============================================================
    fit_sex <- function(sex) {
        
        dat_s <- dat[dat$SEX == sex, ]
        MinLinf <- 150
        MaxLinf <- 250
        
        init <- c(k=0.1, alpha=4, mu=5, sigma=0.5)
        
        fit <- optim(
            init, LL_GL,
            dat=dat_s, MinLinf=MinLinf, MaxLinf=MaxLinf,
            method="L-BFGS-B",
            lower=c(0.01,0.1,0.1,0.01),
            upper=c(1,100,10,10)
        )
        
        Linf <- exp(fit$par[3] + 0.5 * fit$par[4]^2)
        
        if (sex == 1 && hessian_ok(fit$par, dat_s, MinLinf, MaxLinf)) {
            # Female: Hessian SEs
            hess <- optimHess(fit$par, LL_GL,
                              dat=dat_s, MinLinf=MinLinf, MaxLinf=MaxLinf)
            vcov <- solve(hess)
            se_k <- sqrt(vcov[1,1])
            
            grad <- c(Linf, fit$par[4]*Linf)
            se_Linf <- sqrt(grad %*% vcov[3:4,3:4] %*% grad)
            
        } else {
            # Male: Profile SEs
            se_k    <- profile_se_k_fast(fit, dat_s, MinLinf, MaxLinf)
            se_Linf <- profile_se_Linf_fast(fit, dat_s, MinLinf, MaxLinf)
        }
        
        list(
            k = fit$par[1],
            k_SE = se_k,
            Linf = Linf,
            Linf_SE = se_Linf,
            AIC = 2 * fit$value + 2 * length(fit$par)
        )
    }
    
    fem <- fit_sex(1)
    mal <- fit_sex(2)
    
    # ============================================================
    # 6. Final Table
    # ============================================================
    table2 <- data.frame(
        Sex = c("Female","Male"),
        k = c(fem$k, mal$k),
        k_SE = c(fem$k_SE, mal$k_SE),
        Linf = c(fem$Linf, mal$Linf),
        Linf_SE = c(fem$Linf_SE, mal$Linf_SE),
        AIC = c(fem$AIC, mal$AIC)
    )
    
    if (!dir.exists("results/tables"))
        dir.create("results/tables", recursive=TRUE)
    
    write.csv(
        table2,
        "results/tables/table2_gl_model.csv",
        row.names = FALSE
    )
    
    print(table2)
    invisible(table2)
}
test_table2_gl_profile_fast()
