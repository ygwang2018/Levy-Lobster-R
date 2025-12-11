test_table5 <- function(R=200, n_f=50, n_m=50) {
    
    library(dplyr)
    set.seed(123)
    
    ## true values
    kf   <- 0.275;  Lf <- 186.23
    km   <- 0.247;  Lm <- 228.04
    zf   <- 75.055; zm <- 59.947
    gamf <- 12.535; b1f <- -2.087; b2f <- 0.011
    gamm <- 11.33;  b1m <- -1.984; b2m <- 0.010
    
    ## RE parameters
    mlogLf <- log(Lf) - 0.5*0.05^2
    mlogLm <- log(Lm) - 0.5*0.05^2
    shLf <- 200; scLf <- Lf/shLf
    shLm <- 200; scLm <- Lm/shLm
    
    
    ## simulation of one sex
    sim_one <- function(n, sex, k, Linf0, z, gam, b1, b2,
                        type, mlog, sdlog, sh, sc) {
        
        if(n<=0) return(NULL)
        
        if(sex==1) ni <- sample(3:7, n, TRUE) else ni <- sample(2:6, n, TRUE)
        T0 <- sample(10:150, n, TRUE)/365
        
        out <- vector("list", sum(ni))
        id  <- 1L
        
        for(i in seq_len(n)) {
            
            ## Linf draw
            if(type=="fixed") Li <- Linf0 else
                if(type=="lognormal") Li <- rlnorm(1, mlog, sdlog) else
                    Li <- rgamma(1, sh, sc)
                
                PL  <- Li*(1 - exp(-k*T0[i]))
                Tc  <- T0[i]
                
                for(j in seq_len(ni[i])) {
                    
                    a <- (1 - exp(-k*Tc))*(z-1)
                    b <- exp(-k*Tc)*(z-1)
                    if(a<=0 || b<=0) next
                    
                    lam <- rbeta(1,a,b)
                    INCt <- lam*(Li - PL)
                    
                    INCobs <- INCt + rnorm(1,0,1)
                    if(INCobs < 0.1) INCobs <- 0.1
                    
                    mu <- exp(b1 + b2*PL)
                    sh0 <- gam*mu
                    IPt <- rgamma(1, sh0, 1/gam)
                    IPobs <- IPt * rgamma(1,50,1/50)
                    
                    PLobs <- PL + rnorm(1,0,1.5)
                    
                    out[[id]] <- data.frame(
                        LOBSTER=i,
                        SEX=sex,
                        PL=PLobs,
                        INC=INCobs,
                        INT=IPobs*365.25,
                        Tcum=Tc
                    )
                    
                    id <- id+1L
                    PL <- PL + INCt
                    Tc <- Tc + IPt
                }
        }
        bind_rows(out)
    }
    
    
    ## log-likelihood
    loglik <- function(Linf, dat, k, z, gam, b1, b2) {
        
        if(Linf <= max(dat$PL)) return(-1e10)
        LL <- 0
        
        ids <- unique(dat$LOBSTER)
        for(i in ids) {
            d <- dat[dat$LOBSTER==i,]
            Xi <- d$INC; Li <- d$PL
            IPi <- d$INT/365.25
            Tc  <- d$Tcum
            
            for(j in seq_along(Xi)) {
                
                a <- (1-exp(-k*Tc[j]))*(z-1)
                b <- exp(-k*Tc[j])*(z-1)
                if(a<=0 || b<=0) return(-1e10)
                
                den <- Linf - Li[j]
                if(den<=0) return(-1e10)
                
                x <- Xi[j]/den
                if(x<=0 || x>=1) return(-1e10)
                
                logb <- dbeta(x,a,b,log=TRUE) - log(den)
                
                mu <- exp(b1 + b2*Li[j])
                sh0 <- gam*mu
                lg <- dgamma(IPi[j], sh0, scale=1/gam, log=TRUE)
                
                LL <- LL + logb + lg
            }
        }
        LL
    }

  
    ## Linf estimator
    est_Linf <- function(dat, k, z, gam, b1, b2) {
        if(is.null(dat) || nrow(dat)<3) return(NA_real_)
        
        lo <- max(dat$PL) + 1e-3
        up <- lo + 80
        
        out <- suppressWarnings(
            try(optimize(loglik, c(lo,up),
                         dat=dat, k=k, z=z, gam=gam, b1=b1, b2=b2,
                         maximum=TRUE), silent=TRUE)
        )
        
        if(inherits(out,"try-error")) return(NA_real_)
        out$maximum
    }
    
    
    ## summary stats
    ss <- function(est, true) {
        est <- est[is.finite(est)]
        if(length(est)<3) return(c(mean=NA,bias=NA,var=NA,se=NA,rmse=NA))
        
        m  <- mean(est)
        b  <- mean(est - true)
        v  <- var(est)
        s  <- sqrt(v)
        r  <- sqrt(mean((est - true)^2))
        c(mean=m,bias=b,var=v,se=s,rmse=r)
    }
    
    
    ## main simulation
    designs <- c("fixed","lognormal","gamma")
    out <- data.frame()
    
    for(d in designs) {
        
        Fhat <- rep(NA, R)
        Mhat <- rep(NA, R)
        
        for(r in 1:R) {
            
            df <- sim_one(n_f,1,kf,Lf,zf,gamf,b1f,b2f,
                          d, mlogLf,0.05,shLf,scLf)
            
            dm <- sim_one(n_m,2,km,Lm,zm,gamm,b1m,b2m,
                          d, mlogLm,0.05,shLm,scLm)
            
            Fhat[r] <- est_Linf(df, kf, zf, gamf, b1f, b2f)
            Mhat[r] <- est_Linf(dm, km, zm, gamm, b1m, b2m)
        }
        
        sf <- ss(Fhat, Lf)
        sm <- ss(Mhat, Lm)
        
        out <- rbind(out,
                     data.frame(Design=d, Sex="Female", True=Lf,
                                Mean=sf["mean"], Bias=sf["bias"],
                                Variance=sf["var"], SE=sf["se"], RMSE=sf["rmse"]),
                     data.frame(Design=d, Sex="Male", True=Lm,
                                Mean=sm["mean"], Bias=sm["bias"],
                                Variance=sm["var"], SE=sm["se"], RMSE=sm["rmse"])
        )
    }
    
    
    ## save
    if(!dir.exists("results/tables")) dir.create("results/tables", TRUE)
    write.csv(out, "results/tables/table5.csv", row.names=FALSE)
    
    invisible(out)
}
results <- test_table5(R = 200, n_f = 50, n_m = 50)
print(results)
