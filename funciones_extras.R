
##################################################################
# This function computes the true ROC curve based on normal errors
##################################################################

ROC_true_x_new<-function(predict_D,predict_H,sigma_D,sigma_H,p) 
{ 
  if(length(p)==1){
	dif_predict <-predict_H- predict_D  
  	adentro <- qnorm(1-p)*(sigma_H/sigma_D)+dif_predict/sigma_D
  	ROC=(1-pnorm(adentro))
	}
  if(length(p)!=1){
	ROC=rep(NA,length(p))
	for (i  in 1:length(p)){
	
		pe=p[i]
		dif_predict <-predict_H- predict_D  
  		adentro <- qnorm(1-pe)*(sigma_H/sigma_D)+dif_predict/sigma_D
  		ROC[i]=(1-pnorm(adentro))
		}

	}	
  return(ROC)
}		


###################################################################
# This function computes the classical estimator of the ROC curve
# The AUC and YOUDEN INDEX are also returned
###################################################################

CURVA_ROC_CL <- function(residuos_D, residuos_H, predict_D,predict_H,sigma_D,sigma_H,grilla_p){
	residuos_D_std<- residuos_D/sigma_D
  	residuos_H_std<- residuos_H/sigma_H

  	ROC_CL_grilla<- rep(NA,length= length(grilla_p)) 
 
	dif_predict <-predict_H- predict_D 
	for(i in 1:length(grilla_p)){
			p<-grilla_p[i]
      		cuantil_H<- unname(quantile(residuos_H_std,probs=1-p))
			punto<- cuantil_H*(sigma_H/sigma_D)+dif_predict/sigma_D 
			punto <- as.numeric(punto)
			ROC_CL_grilla[i]<- 1-mean(1*(residuos_D_std<=punto))
    	}
    
	AUC_CL<-mean(ROC_CL_grilla)
	YOUDEN_CL <- max(ROC_CL_grilla- grilla_p)
	return(list(ROC_CL_grilla=ROC_CL_grilla, AUC_CL=AUC_CL, YOUDEN_CL=YOUDEN_CL))
}


###################################################################
# This function computes the robust estimator of the ROC curve
# The AUC and YOUDEN INDEX are also returned
###################################################################

CURVA_ROC_ROB<-function(residuos_D, residuos_H, predict_D,predict_H,sigma_D,sigma_H, eta=2.5, grilla_p){
	residuos_D_std<- residuos_D/sigma_D
  	residuos_H_std<- residuos_H/sigma_H

  	ROC_ROB_grilla<- rep(NA,length= length(grilla_p)) 

  	prop_D <-proporcion(residuos_D_std,eta=eta)
  	prop_H<- proporcion(residuos_H_std,eta=eta)
	
	## Calculo pesos
 
  	if(prop_D$t_ene==0){
  		pesos_D=rep(1,length(residuos_D))
  	}
  	if(prop_D$t_ene!=0){
  		pesos_D<- whard(abs(residuos_D_std)/prop_D$t_ene)
  	}

  	if(prop_H$t_ene==0){
  		pesos_H=rep(1,length(residuos_H))
  	}
  
	if(prop_H$t_ene!=0){
  		pesos_H<- whard(abs(residuos_H_std)/prop_H$t_ene)
  	}
	 

      dif_predict <-predict_H-predict_D 

      for(i in 1:length(grilla_p)){
           	pe=grilla_p[i]
           	cuantil_H <- unname(weighted.fractile(residuos_H_std,pesos_H,1-pe))
           	punto<- cuantil_H*(sigma_H/sigma_D)+dif_predict/sigma_D
		punto <- as.numeric(punto)
		ROC_ROB_grilla[i]<- 1-fdaw(residuos_D_std,pesos_D,punto) 
          }
      
	AUC_ROB<- mean(ROC_ROB_grilla)
	YOUDEN_ROB <- max(ROC_ROB_grilla- grilla_p)
      	return(list(ROC_ROB_grilla=ROC_ROB_grilla, AUC_ROB=AUC_ROB, YOUDEN_ROB=YOUDEN_ROB))
}

#################################################################
# These functions are used internally by CURVA_ROC_ROB
#################################################################

#################################################################
# WEIGHTED EMPIRICAL DISTRIBUTION
##############################################################


fdaw <-function(x,wei,t){
  wei=wei/sum(wei)
	fdaw<- sum(wei[x <= t])
	fdaw	
}

#################################################################
# WEIGHTED QUANTILE
##############################################################


weighted.fractile<-function (y, w, p){
    w <- w/sum(w)
    a <- 1 - p
    b <- p
    ox <- order(y)
    y <- y[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0, w))
    up <- sum(w) - low
    df <- a * low - b * up
    repeat {
        if (df[k] < 0) 
            k <- k + 1
        else if (df[k] == 0) 
            return((w[k] * y[k] + w[k - 1] * y[k - 1])/(w[k] + 
                w[k - 1]))
        else return(y[k - 1])
    }
}


#################################################################
# These functions are used internally by CURVA_ROC_ROB
# To compute the cut-off values
#################################################################

#################################################################
Phimas<-function(x){
	Phimas<- 1-2*(1-pnorm(abs(x)))
	Phimas
}

######################################################################
parte_plus<-function(x){ 
   parte_plus<- max(x,0)
   parte_plus}

#################################################################
proporcion<-function(res,eta=2.5){
  ene<- length(res)
  abs.res<- abs(res)
	abs.resor<- sort(abs.res)
  indica_eta<- 1*(abs.resor< eta)
  i0<- sum(indica_eta)
  
 
  if (i0 == ene){
  	proporcion<- list(i0=ene,d_ene=0,i_ene=ene,t_ene=0)
  	return(proporcion)}
  
  indices<- seq(1,ene,1)*(1-indica_eta)
  diferencias<- as.matrix(Phimas(abs.resor[indica_eta==0])-((indices[indica_eta==0]-1)/ene),nrow=ene,ncol=1)
  diferencias_plus<- apply(diferencias,1,parte_plus)
  
 
  i_maximo <- which.max(diferencias_plus)
  d_ene <- diferencias_plus[i_maximo]
  
  i_ene <-  ene-floor(ene*d_ene)
  t_ene<-   abs.resor[i_ene]
  
  proporcion<- list(i0=i0,d_ene=d_ene,i_ene=i_ene,t_ene=t_ene)
  return(proporcion)
}

#################################################################
# Weight function w to be used to downweight large residuals
##############################################################

whard<-function(x){
	whard<-ifelse(abs(x)<1,1,0)
	whard	
}


tukey.weight <- function(a) {
  tmp <- 6*(1-a^2)^2
  tmp[ abs(a) > 1 ] <- 0
  return(tmp)
}
