#########################################################
# Libraries and codes needed
#########################################################

library(fda)                     # splines
library(robustbase)              # lmrob

source("FLM_fit.R")
source("generar_ROC.R")
source("funciones_extras.R")

######################################################
# FUNCTION USED TO PERFORM THE SIMULATIONS
######################################################
     

simulacion_ROC<- function(Nrep=1000, nn = c(300,300), grillat= seq(0, 1, length = 100), modelo =1, 
                   conta=list(type = 'C0',ratio=0, media_xx = NA , ee_c_mu =NA),
		   slope=c(1.5,1), alfas = c(2,1.5), ee_sd = c(2,1),
                   criterio_seleccion = "bic1", funcionx0 = "cero",rango_D, rango_H) {


lt <- length(grillat)
media_ee_c <- conta$ee_c_mu[1]
media_xx <- conta$xx_mu[1]

################################################################
# We first compute the true ROC
################################################################
ttt   <- grillat
phi_1   <-  sqrt(2) * sin(0.5*pi * ttt) 

dt <- ttt[2]-ttt[1]

beta_D <- slope[1]* beta_FUN(ttt)
beta_H <- slope[2]* beta_FUN(ttt)

coef_beta_D <- slope[1] * c(1,1)
coef_beta_H <- slope[2] * c(1,1)

X0t <- rep(0, length(ttt))

    if (funcionx0=="seno"){
           X0t <- x0(ttt)
	}

#########################
# Computes <X0, beta_D>
#########################

corro_D <- as.numeric(t(X0t)%*% beta_D*dt )
corro_H <- as.numeric( t(X0t)%*% beta_H*dt)

###############################################################
# DEFINE THE POSSIBLE VALUES WHERE TO EVALUATE THE ROC_X
# AS  Xz= X0t+ Z*phi_1 
# with Z varying between -2*sqrt(lambda1) and 2*sqrt(lambda1)
# lambda1 being the first eigenvalue of the covariance operator
################################################################
lambda1 <- 4
zetas <- seq(-2*sqrt(lambda1), 2*sqrt(lambda1),  by=0.05)
lzetas  <-  length(zetas)

grilla_p <- seq(0.01, 0.99,  by=0.01)
lpes <-  length(grilla_p)

#############################################
# DEFINE Xz
############################################
phi_1_t   <-  sqrt(2) * sin(0.5*pi * grillat)

Xz <- matrix(NA, lzetas*lt,nrow=lzetas,ncol=lt) 

for(k in 1:lzetas){
	Xz[k,] <-X0t + zetas[k] * phi_1_t 
}

#############################################
# PLOT Xz
############################################

par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
matplot(grillat, t(Xz), xlab="t", ylab="", type = 'l', col = "grey", lty = 1)

##############################################
# COMPUTE THE TRUE ROC_X
#############################################
	 
ROC_normal_x<- matrix(rep(NA,lzetas*lpes),nrow=lzetas,ncol=lpes )
YOUDEN_normal <- AUC_normal <- rep(NA, lzetas)

for(k in 1:lzetas){
	# PORQUE si X_Z= X0+ z *phi_1 entonces
	# < X_z, beta_D > = corro_D + z* slope[1]

      predict_D <- alfas[1]+ corro_D +  zetas[k] * slope[1]

      predict_H <- alfas[2]+ corro_H +  zetas[k] * slope[2]
 
	sigma_D <- ee_sd[1]
	sigma_H <- ee_sd[2]

	ROC_normal_x[k,]<- ROC_true_x_new(predict_D,predict_H,sigma_D,sigma_H, grilla_p) 
        AUC_normal[k]<- mean(ROC_normal_x[k,])
	YOUDEN_normal[k] <- max(ROC_normal_x[k,]- grilla_p)
 }

#########################################
# CHANGE THE FOLDER TO SAVE THE RESULTS
#########################################
carpeta <-paste("SALIDAS_n",ene,"_mod_",modelo,"_funcionx0_", funcionx0, sep="")    

if (!dir.exists(carpeta)) dir.create(carpeta)

setwd(carpeta)

#####################################################
# NAMES OF THE FILES
######################################################
 

ARCHIVO_D_ROB <- paste("ALFA_BETA_ROB_D_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")   
ARCHIVO_H_ROB <- paste("ALFA_BETA_ROB_H_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")   
ARCHIVO_D_CL <- paste("ALFA_BETA_CL_D_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")   
ARCHIVO_H_CL <- paste("ALFA_BETA_CL_H_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")   
if(conta$type== 'C0'){
	ARCHIVO_D_ROB <-paste("ALFA_BETA_ROB_D_n",ene,"_mod_",modelo,"_cont_",conta$type,".txt",sep="")
	ARCHIVO_H_ROB <-paste("ALFA_BETA_ROB_H_n",ene,"_mod_",modelo,"_cont_",conta$type,".txt",sep="")   
	ARCHIVO_D_CL  <-paste("ALFA_BETA_CL_D_n",ene,"_mod_",modelo,"_cont_",conta$type,".txt",sep="")
	ARCHIVO_H_CL  <-paste("ALFA_BETA_CL_H_n",ene,"_mod_",modelo,"_cont_",conta$type,".txt",sep="")   
  }
 
###################################################
# FILES WITH SUMMARY MEASURES FOR AUC y ROC
###################################################

ARCHIVOAUC_ROB <-paste("AUC_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")
ARCHIVOMSE_ROB <-paste("MSE_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")
ARCHIVOKS_ROB <-paste("KS_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")


ARCHIVOAUC_CL <-paste("AUC_CL_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")
ARCHIVOMSE_CL <-paste("MSE_CL_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")
ARCHIVOKS_CL <-paste("KS_CL_n",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".txt",sep="")


###################################################
# SHORTEN THE NAME FOR CLEAN SAMPLES
###################################################

if(conta$type== 'C0'){
	ARCHIVOAUC_ROB <-paste("AUC_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type,".txt",sep="")
	ARCHIVOMSE_ROB <-paste("MSE_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type, ".txt",sep="")
	ARCHIVOKS_ROB <-paste("KS_ROB_n",ene,"_mod_",modelo,"_cont_",conta$type, ".txt",sep="")


	ARCHIVOAUC_CL <-paste("AUC_CL_n",ene,"_mod_",modelo,"_cont_",conta$type, ".txt",sep="")
	ARCHIVOMSE_CL <-paste("MSE_CL_n",ene,"_mod_",modelo,"_cont_",conta$type, ".txt",sep="")
	ARCHIVOKS_CL <-paste("KS_CL_n",ene,"_mod_",modelo,"_cont_",conta$type, ".txt",sep="")
}

########################################################################################
# SAVE IN ARRAYS TO SAVE THE SIMULATION RESULTS IN AN Rdata
########################################################################################


alfa_D_CL=rep(NA,Nrep)
alfa_D_ROB=rep(NA,Nrep)
beta_D_CL=matrix(NA, Nrep,lt)
beta_D_ROB=matrix(NA, Nrep,lt)

sigma_D_CL=rep(NA,Nrep)
sigma_D_ROB=rep(NA,Nrep)

alfa_H_CL=rep(NA,Nrep)
alfa_H_ROB=rep(NA,Nrep)
beta_H_CL=matrix(NA, Nrep,lt)
beta_H_ROB=matrix(NA, Nrep,lt)

sigma_H_CL=rep(NA,Nrep)
sigma_H_ROB=rep(NA,Nrep)

iteracion =rep(NA,Nrep)

orden_D_CL=rep(NA,Nrep)
orden_D_ROB=rep(NA,Nrep)

orden_H_CL=rep(NA,Nrep)
orden_H_ROB=rep(NA,Nrep)

####################################
# SAVE THE ROC ESTIMATES IN MATRICES
####################################
dimension_Arreglo_est <- c(Nrep,lzetas,lpes)

ROC_CL_est <- array(rep(NA,Nrep*lzetas*lpes),dim= dimension_Arreglo_est)
YOUDEN_CL_est<- AUC_CL_est <- array(rep(NA,Nrep*lzetas),dim= dimension_Arreglo_est[1:2])

ROC_ROB_est <- array(rep(NA,Nrep*lzetas*lpes),dim= dimension_Arreglo_est)
YOUDEN_ROB_est<- AUC_ROB_est <- array(rep(NA,Nrep*lzetas),dim= dimension_Arreglo_est[1:2])

KS_YOUDEN_CL_est  <-KS_AUC_CL_est  <-  KS_ROC_CL_est <- rep(NA,Nrep)
MSE_YOUDEN_CL_est  <- MSE_AUC_CL_est <- MSE_ROC_CL_est <- rep(NA,Nrep)


KS_YOUDEN_ROB_est  <-KS_AUC_ROB_est  <-  KS_ROC_ROB_est <- rep(NA,Nrep)
MSE_YOUDEN_ROB_est  <- MSE_AUC_ROB_est <- MSE_ROC_ROB_est <- rep(NA,Nrep)
  
 
#################################################
# START THE REPLICATIONS
##################################################
t1=Sys.time()

##################################################
## SOME SAMPLES MAY BE DISCARDED...
##################################################

done <- 0
desde<- 1
hasta <- Nrep
iter <- desde


while(done <= (hasta - desde)){
	print(iter)
	#########################################
	# GENERATE THE SAMPLES
	###########################################

 	smpl <- generar(iter,funcionx0=funcionx0, nn=nn, modelo=modelo, grillat=grillat, alfas=alfas,  slope=slope, ee_sd=ee_sd, cont=conta)
	## Integration step
	dt <-  grillat[2]-grillat[1]

	t_D = t_H = grillat

	y_D=smpl$y_D
	x_D=smpl$x_D
        out_D= smpl$o_D

	y_H=smpl$y_H
	x_H=smpl$x_H
        out_H= smpl$o_H

	#########################################
	# ROBUST REGRESSION ESTIMATORS
	###########################################
 
	aa_D= try(FLMBsplines(y_D, x_D, t_D, range_freq = rango_D,
                         norder = 4,
                         fLoss = 'lmrob', criterion = criterio_seleccion,
                         trace = FALSE))
	if (!inherits(aa_D, 'try-error')) {
		aa_H= try(FLMBsplines(y_H, x_H, t_H, range_freq = rango_H,
                         norder = 4,
                         fLoss = 'lmrob', criterion = criterio_seleccion,
                         trace = FALSE))

		if (!inherits(aa_H, 'try-error')) {

			iter1 = done+1

			alfa_D_ROB[iter1]  = aa_D$fit$intercept
			beta_D_ROB[iter1,] = aa_D$fit$slope_fun
			sigma_D_ROB[iter1] = aa_D$fit$scale
			orden_D_ROB[iter1] = aa_D$freq

			alfa_H_ROB[iter1]  = aa_H$fit$intercept
			beta_H_ROB[iter1,] = aa_H$fit$slope_fun
			sigma_H_ROB[iter1] = aa_H$fit$scale
			orden_H_ROB[iter1] = aa_H$freq

			vector_D_ROB<- c(iter, aa_D$freq, aa_D$fit$scale, aa_D$fit$intercept, aa_D$fit$slope_fun)
			vector_H_ROB<- c(iter, aa_H$freq, aa_H$fit$scale, aa_H$fit$intercept, aa_H$fit$slope_fun)
	
			lvec_D=length(vector_D_ROB)
 			lvec_H=length(vector_H_ROB)
   
			write(t(vector_D_ROB),file=ARCHIVO_D_ROB,ncolumns=lvec_D,append=T)
			write(t(vector_H_ROB),file=ARCHIVO_H_ROB,ncolumns=lvec_H,append=T)

			predicho_DATA_D_ROB= aa_D$fit$intercept+ x_D%*% aa_D$fit$slope_fun*dt
			predicho_DATA_H_ROB= aa_H$fit$intercept+ x_H%*% aa_H$fit$slope_fun*dt

			residuo_D_ROB= y_D- predicho_DATA_D_ROB
			residuo_H_ROB= y_H- predicho_DATA_H_ROB

			###########################################
			# COMPUTES THE ROBUST ESTIMATOR OF THE ROC
			###########################################
 			for(k in 1:lzetas){

				X_z <- Xz[k,]
      			        predicho_D_ROB <- aa_D$fit$intercept+  +  X_z%*% aa_D$fit$slope_fun*dt
     				predicho_H_ROB <- aa_H$fit$intercept+  +  X_z%*% aa_H$fit$slope_fun*dt
				
				calculo_ROB <- CURVA_ROC_ROB(residuos_D=residuo_D_ROB, residuos_H= residuo_H_ROB, 
							predict_D=predicho_D_ROB,predict_H=predicho_H_ROB,
							sigma_D=aa_D$fit$scale, sigma_H=aa_H$fit$scale, eta=2.5, grilla_p)
				ROC_ROB_est[iter1,k,] <- calculo_ROB$ROC_ROB_grilla
  				AUC_ROB_est[iter1,k] <- calculo_ROB$AUC_ROB
				YOUDEN_ROB_est[iter1,k] <- calculo_ROB$YOUDEN_ROB
			}
  
			#########################################
			# SUMMARY MEASURES
			###########################################
 
			KS_AUC_ROB_est[iter1] <- max(abs(AUC_normal- AUC_ROB_est[iter1,]))
 			KS_ROC_ROB_est[iter1] <- max(abs(ROC_normal_x- ROC_ROB_est[iter1,,]))
			KS_YOUDEN_ROB_est[iter1] <- max(abs(YOUDEN_normal- YOUDEN_ROB_est[iter1,]))
			MSE_AUC_ROB_est[iter1] <- mean((AUC_normal- AUC_ROB_est[iter1,])^2)
 			MSE_ROC_ROB_est[iter1] <- mean((ROC_normal_x- ROC_ROB_est[iter1,,])^2)
			MSE_YOUDEN_ROB_est[iter1] <- mean((YOUDEN_normal- YOUDEN_ROB_est[iter1,])^2)


			vector_ROB<- c(iter, aa_D$freq,aa_H$freq, KS_ROC_ROB_est[iter1], KS_AUC_ROB_est[iter1],KS_YOUDEN_ROB_est[iter1])
 			lvec <- length(vector_ROB)
			write(t(vector_ROB),file=ARCHIVOKS_ROB,ncolumns=lvec,append=T) 

			vector_ROB<- c(iter, aa_D$freq,aa_H$freq, MSE_ROC_ROB_est[iter1], MSE_AUC_ROB_est[iter1],MSE_YOUDEN_ROB_est[iter1])
 			lvec <- length(vector_ROB)
			write(t(vector_ROB),file=ARCHIVOMSE_ROB,ncolumns=lvec,append=T) 

			vector_ROB<- c(iter, aa_D$freq,aa_H$freq, AUC_ROB_est[iter1,])
 			lvec <- length(vector_ROB)
			write(t(vector_ROB),file=ARCHIVOAUC_ROB,ncolumns=lvec,append=T) 


			#########################################
			# CLASSICAL REGRESSION ESTIMATORS 
			###########################################
 
			bb_D=FLMBsplines(y_D, x_D, t_D, range_freq = rango_D,
                         norder = 4,
                         fLoss = 'ls', criterion =criterio_seleccion,
                         trace = FALSE)

			bb_H=FLMBsplines(y_H, x_H, t_H, range_freq = rango_H,
                         norder = 4,
                         fLoss = 'ls', criterion =criterio_seleccion,
                         trace = FALSE)

		
		 	alfa_D_CL[iter1]  = bb_D$fit$intercept
			beta_D_CL[iter1,] = bb_D$fit$slope_fun
			sigma_D_CL[iter1] = bb_D$fit$scale
			orden_D_CL[iter1] = bb_D$freq

			alfa_H_CL[iter1]  = bb_H$fit$intercept
			beta_H_CL[iter1,] = bb_H$fit$slope_fun
			sigma_H_CL[iter1] = bb_H$fit$scale
			orden_H_CL[iter1] = bb_H$freq

			vector_D_CL<- c(iter, bb_D$freq, bb_D$fit$scale, bb_D$fit$intercept, bb_D$fit$slope_fun)
			vector_H_CL<- c(iter, bb_H$freq, bb_H$fit$scale, bb_H$fit$intercept, bb_H$fit$slope_fun)
	
			lvec_D=length(vector_D_CL)
 			lvec_H=length(vector_H_CL)
   
			write(t(vector_D_CL),file=ARCHIVO_D_CL,ncolumns=lvec_D,append=T)
			write(t(vector_H_CL),file=ARCHIVO_H_CL,ncolumns=lvec_H,append=T)

			predicho_DATA_D_CL= bb_D$fit$intercept+ x_D%*% bb_D$fit$slope_fun*dt
			predicho_DATA_H_CL= bb_H$fit$intercept+ x_H%*% bb_H$fit$slope_fun*dt

			residuo_D_CL= y_D- predicho_DATA_D_CL
			residuo_H_CL= y_H- predicho_DATA_H_CL

 			###############################################
			# COMPUTES THE CLASSICAL ESTIMATOR OF THE ROC
			###############################################
			for(k in 1:lzetas){
 
				X_z <- Xz[k,]
     				predicho_D_CL <- bb_D$fit$intercept+  +  X_z%*% bb_D$fit$slope_fun*dt
     				predicho_H_CL <- bb_H$fit$intercept+  +  X_z%*% bb_H$fit$slope_fun*dt
				
				calculo_CL <- CURVA_ROC_CL(residuos_D=residuo_D_CL, residuos_H= residuo_H_CL, 
							predict_D=predicho_D_CL,predict_H=predicho_H_CL,
							sigma_D=bb_D$fit$scale, sigma_H=bb_H$fit$scale,  grilla_p)
				ROC_CL_est[iter1,k,] <- calculo_CL$ROC_CL_grilla
  				AUC_CL_est[iter1,k] <- calculo_CL$AUC_CL
				YOUDEN_CL_est[iter1,k] <- calculo_CL$YOUDEN_CL
			}
  
			#########################################
			# SUMMARY MEASURES
			###########################################
 
			KS_AUC_CL_est[iter1] <- max(abs(AUC_normal- AUC_CL_est[iter1,]))
 			KS_ROC_CL_est[iter1] <- max(abs(ROC_normal_x- ROC_CL_est[iter1,,]))
			KS_YOUDEN_CL_est[iter1] <- max(abs(YOUDEN_normal- YOUDEN_CL_est[iter1,]))
			MSE_AUC_CL_est[iter1] <- mean((AUC_normal- AUC_CL_est[iter1,])^2)
 			MSE_ROC_CL_est[iter1] <- mean((ROC_normal_x- ROC_CL_est[iter1,,])^2)
			MSE_YOUDEN_CL_est[iter1] <- mean((YOUDEN_normal- YOUDEN_CL_est[iter1,])^2)

 
			#########################################
			# SAVE THE RESULTS
			###########################################
 
			vector_CL <- c(iter, bb_D$freq, bb_H$freq, KS_ROC_CL_est[iter1], KS_AUC_CL_est[iter1],KS_YOUDEN_CL_est[iter1])
 			lvec <- length(vector_CL)
			write(t(vector_CL),file=ARCHIVOKS_CL,ncolumns=lvec,append=T) 

			vector_CL<- c(iter, bb_D$freq, bb_H$freq, MSE_ROC_CL_est[iter1], MSE_AUC_CL_est[iter1],MSE_YOUDEN_CL_est[iter1])
 			lvec <- length(vector_CL)
			write(t(vector_CL),file=ARCHIVOMSE_CL,ncolumns=lvec,append=T) 

			vector_CL<- c(iter, bb_D$freq, bb_H$freq, AUC_CL_est[iter1,])
 			lvec <- length(vector_CL)
			write(t(vector_CL),file=ARCHIVOAUC_CL,ncolumns=lvec,append=T) 


			iteracion[iter1] = iter

			## Update number of successful iterations
			done <- done + 1
		}
	}

	print(c(iter, done))
	iter <- iter + 1
}

 t2=Sys.time()

tiempo=t2-t1

tiempo
archivo <- paste("simu-ROC-CL-ROB-n-",ene,"_mod_",modelo,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  media_ee_c, "_xxmu_", media_xx,".RData", sep="")
   
   
if(conta$type== 'C0'){
archivo <- paste("simu-ROC-CL-ROB-n-",ene,"_mod_",modelo,"_cont_",conta$type,".RData", sep="")
}


}


