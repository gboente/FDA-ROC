############################################################################
# Estimators for FLM as defined in Boente et al. (2020)
# The code is a slight modification of the one given 
# by Matias Salibian Barrera for Functional Partial Linear Models
# at https://github.com/msalibian/RobustFPLM
############################################################################

FLMBsplines <- function(y, x, t, range_freq = range_default,
                         norder = 4,
                         fLoss = "lmrob", criterion = "bic1",
                         trace = FALSE) {

    ## Some Setup
    opt <-   freq_opt <- fit_opt <- Inf
    n <- length(y)
    range_default <- floor(max(n^(1 / 5), norder)):
        floor(2 * (norder + n^(1 / 5)))

    ## Loop
         for (freq in range_freq) {
            fit <- FLMBsplines_fit(y, x,   t, freq,   norder, fLoss)
            val <- fit$value
            scl <- fit$scale
            crt <- goodness1(n, scl, val, freq,  criterion)
            if (crt < opt) {
                opt <- crt
                freq_opt <- freq
                fit_opt <- fit
            }
            if (trace) print(c("freq" = freq, "crit" = crt))
        }
 

   
    dt <- min(diff(t))
    fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt + fit_opt$intercept)

    return(list(fit = fit_opt, freq = freq_opt))
}


FLMBsplines_fit <- function(y, x, t, freq,  norder, fLoss) {

    ## Integration step
    dt <- min(diff(t)) # width of grid
    xcenter <- x

    ## SPLINE BASIS 

    grilla_spl <- t # seq(0, 1, length = length(t))
    nodos_spl <- seq(min(grilla_spl), max(grilla_spl),
                     length = freq - norder + 2)
    base_spl <- create.bspline.basis(
        rangeval = range(t),
        norder = norder,
        breaks = nodos_spl
    )
    beta_spl <- getbasismatrix(grilla_spl, base_spl)
    cov_dec <- beta_spl[, 1:freq]

    ## Estimated coefficients (by row)
    xx_coef <- xcenter %*% cov_dec * dt

    ## Parameter estimation
    est <- minimizo(y, xx_coef,   fLoss  )

    est$slope_fun <- cov_dec %*% est$slope
    est$intercept <- est$alfa
    return(est)
}

##################################################################
# This function is used internally by FLMBsplines_fit 
##################################################################

minimizo <- function(y, x_coef,   fLoss, 
                     rob.control = lmrob.control(
                         trace.level = 0,
                         nResample = 5000, 
                         tuning.psi = 3.443689, # 85% eff
                         subsampling = "simple",
                         rel.tol = 1e-5, 
                         refine.tol = 1e-5,
                         k.max = 2e3,
                         maxit.scale = 2e3,
                         max.it = 2e3
                     )) {

    ## Design matrix
    X <- x_coef 

    if (!(fLoss %in% c("ls", "huang", "lmrob"))) {
        stop("Invalid fLoss. Should be one of \'ls\' or \'huang\' or \'lmrob\' ")
    }

    ## Minimization menu
    switch(fLoss,
           ls = {
               fit <- lm(y ~ X )
               cf <- fit$coef
               vv <- sum(fit$res^2) # 1
               ss <- sqrt(mean(fit$res^2))
           },
           lmrob = {
               fit <- lmrob(y ~ X, control = rob.control)
               if (fit$init.S$converged) {
                   cf <- fit$coef
                   ss <- fit$scale
                   vv <- sum(Mpsi(fit$res / ss,
                                  cc = rob.control$tuning.psi,
                                  psi = rob.control$psi, deriv = -1
                                  ))
               } else {
                   stop("S-estimator did not converge.")
               }
           },
           stop("Invalid fLoss.")
           )

    spl_intercept <- cf[1]
    slope_par <- cf[-1]

    return(list(
        alfa = spl_intercept,
        slope = slope_par,
        value = vv,
        scale = ss, 
       residuo=fit$resid
    ))
}



##################################################################
# This function is used internally by FLMBsplines 
# It computes the goodness criterion to be minimized
##################################################################

goodness1 <- function(nn, scl, val,  freq, criterion) {
  switch(criterion,
    aic = log(scl^2 * val / nn)+ 2 * freq / nn,
    bic = log(scl^2 * val / nn) + freq * log(nn) / (2 * nn),
    bic1 = log(scl^2 * val / nn) +   freq  * log(nn) / (nn) 
  )
}