## ------------------------------------------------------
## Generates clean and contaminated samples 
## ------------------------------------------------------
## INPUT
##
## iter  : seed
## nn    : sample size for each population
## modelo: Two possible models according to Inacio et al (2012)
## lt    : grid size (functional term) 
## alfas : intercept for each population
## alfas[1] : intercept for diseased
## alfas[2] : intercept for healthy
## slope : constant multiplying the slope beta for each population
## slope[1] : constant for diseased
## slope[2] : constant for healthy
## ee_sd : error sd for each population
## ee_sd[1] : error sd for diseased
## ee_sd[2] : error sd for healthy
## cont  : list of contamination parameters
##           type  : 'C0','C1H','C2H','C1D','C2D','C1HD','C2HD',
##           ratio : contamination ratio, the same for both if both are contaminated
##           xx_col: contaminated column, the same for both if both are contaminated
##           xx_mu : mean, [1] diseased [2] healthy when contaminating
##           xx_sd : sd, [1] diseased [2] healthy when contaminating
##           ee_c_mu : mean vector [1] diseased [2] healthy, corresponds to a shift when contaminating
##           ee_c_sd : sd when contaminating
## ------------------------------------------------------
## OUTPUT
##
## y    : response
## x    : functional observations
## xi   : xx's Fourier coefficients
## xic  : xc's Fourier coefficients (after contamination)
## e    : errors
## phi  : true L^2 basis
## o_D    : outlier indices for diseased
## o_H    : outlier indices for healthy
## ------------------------------------------------------
x0<-function(t){t+sin(t)}

beta_FUN <- function(t){ sqrt(2)*(sin(pi*t/2)+ sin(1.5*pi*t))}

generar <- function (iter, funcionx0="cero", nn = c(300,300), modelo=1, grillat= seq(0, 1, length = 100), alfas=c(2,1.5),  slope=c(1.5,1),
                     ee_sd=c(2,1), cont) {

    ## Set seed for reproducibility
    set.seed(iter)

    nn_D <- nn[1]
    nn_H <- nn[2]

    ## Define the orthonormal L^2 Basis (by columns)
    ttt   <- grillat
    phi   <- function (j, t) {sqrt(2) * sin(j*0.5* pi * t )}
    phi_1   <-  phi(1,ttt) 
    phi_2   <-  phi(3,ttt) 
    phi_3   <-  phi(5,ttt) 
    phi_4   <-  phi(7,ttt) 
    phi_5   <-  phi(9,ttt) 
    phi_j   <- cbind(phi_1, phi_2, phi_3,phi_4,phi_5)
    autovalores <- c(4, 1, 0.5,0.1,0.05)

    dt <- ttt[2]-ttt[1]
     
    beta_D <- slope[1]* beta_FUN(ttt)
    beta_H <- slope[2]* beta_FUN(ttt)

    coef_beta_D <- slope[1] * c(1,1,0,0,0)
    coef_beta_H <- slope[2] * c(1,1,0,0,0)

    ## XX fourier coefficients (by row)
    xx_D_coef <- xc_D_coef <- t(replicate(nn_D, rnorm(5, 0, sqrt(autovalores))))
    xx_H_coef <- xc_H_coef <- t(replicate(nn_H, rnorm(5, 0, sqrt(autovalores))))


    ## Contamination 
    out_D <- out_H <- NULL
    ee_D  <- ee_sd[1]*rnorm(nn_D)
    ee_H  <- ee_sd[2]*rnorm(nn_H)

    ## ...in the response Healthy
    if (cont$type == 'C1H') {
        out_H <- rbinom(nn_H, 1, cont$ratio)
        ee_H[out_H == 1] <- cont$ee_c_mu[2]+  cont$ee_c_sd[2]* rnorm(sum(out_H))
    }
    ## ...in the response DISEASED

    if (cont$type == 'C1D') {
        out_D <- rbinom(nn_D, 1, cont$ratio)
        ee_D[out_D == 1] <- cont$ee_c_mu[1]+  cont$ee_c_sd[1]* rnorm(sum(out_D))

    }

    if (cont$type == 'C1HD') {
        out_H <- rbinom(nn_H, 1, cont$ratio)
        ee_H[out_H == 1] <- cont$ee_c_mu[2]+  cont$ee_c_sd[2]* rnorm(sum(out_H))
        out_D <- rbinom(nn_D, 1, cont$ratio)
        ee_D[out_D == 1] <- cont$ee_c_mu[1]+  cont$ee_c_sd[1]* rnorm(sum(out_D))

    }

    ## ...both in the response and carriers for healthy
    if (cont$type == 'C2H') {
        out_H <- rbinom(nn_H, 1, cont$ratio)
        ee_H[out_H == 1] <- cont$ee_c_mu[2]+  cont$ee_c_sd[2]* rnorm(sum(out_H))
        xc_H_coef[out_H == 1, cont$xx_col] <- cont$xx_mu[2] + cont$xx_sd  * rnorm(sum(out_H))

    }
     if (cont$type == 'C2D') {
        out_D <- rbinom(nn_D, 1, cont$ratio)
        ee_D[out_D == 1] <- cont$ee_c_mu[1]+  cont$ee_c_sd[1]* rnorm(sum(out_D))
        xc_D_coef[out_D == 1, cont$xx_col] <- cont$xx_mu[1] + cont$xx_sd  * rnorm(sum(out_D))

    }
    if (cont$type == 'C2HD') {
        out_H <- rbinom(nn_H, 1, cont$ratio)
        ee_H[out_H == 1] <- cont$ee_c_mu[2]+  cont$ee_c_sd[2]* rnorm(sum(out_H))
        xc_H_coef[out_H == 1, cont$xx_col] <- cont$xx_mu[2] + cont$xx_sd[2] * rnorm(sum(out_H))
        out_D <- rbinom(nn_D, 1, cont$ratio)
        ee_D[out_D == 1] <- cont$ee_c_mu[1]+  cont$ee_c_sd[1]* rnorm(sum(out_D))
        xc_D_coef[out_D == 1, cont$xx_col] <- cont$xx_mu[1] + cont$xx_sd[1] * rnorm(sum(out_D))

    }


    ## Functional data (by row)
    X0t <- rep(0, length(ttt))

    if (funcionx0=='seno'){
           X0t <- x0(ttt)
	}
    xx_H_fun <- X0t+ xc_H_coef %*% t(phi_j)  
    xx_D_fun <- X0t+ xc_D_coef %*% t(phi_j)  
     
    corro_D <- as.numeric(t(X0t)%*% beta_D*dt) # CALCULA <X0, beta_D>
    corro_H <- as.numeric( t(X0t)%*% beta_H*dt )

    regre_D <- alfas[1]+ corro_D  + as.vector(xc_D_coef %*% coef_beta_D )
    regre_H <- alfas[2]+ corro_H + as.vector(xc_H_coef %*% coef_beta_H)

    if(modelo==2){
	xx_D_fun <- xx_D_fun^2
	regre_D <- alfas[1]+ as.numeric( xx_D_fun%*% beta_D*dt) 
    }

    yy_D  <- regre_D  + ee_D
    yy_H  <- regre_H   + ee_H
    
    
    return(list(y_D   = yy_D, y_H   = yy_H,
                x_D   = xx_D_fun, x_H   = xx_H_fun, 
                beta_D  =  beta_D, beta_H  =  beta_H,
                fx_D  = regre_D ,
		    fx_H  = regre_H ,
                o_D   = out_D, o_H   = out_H))

}


