rm(list = ls())

library(fda)                     # splines
library(robustbase)              # lmrob

source("FLM_fit.R")
source("generar_ROC.R")
source("funciones_extras.R")
source("function_simulation_ROC.R")
######################################################
# PARAMETROS SIMULACION
######################################################

Nrep=1000

ene <- 300
nn = c(ene,ene)


nn_D <- nn[1]
nn_H <- nn[2]


lt = 1e2
grillat= seq(0, 1, length = lt)

modelo =1
slope=c(1.5,1)
alfas = c(2,1.5)
ee_sd = c(2,1)

criterio_seleccion = "bic1"
funcionx0 = "cero"


menorp_D = max(floor(nn_D^(1/5)/2), 4)
mayorp_D =  floor(8+2*nn_D^(1/5))


menorp_H = max(floor(nn_H^(1/5)/2), 4)
mayorp_H =  floor(8+2*nn_H^(1/5))

rango_D=c(menorp_D:mayorp_D)
rango_H=c(menorp_H:mayorp_H)

##################################
# CLEAN SAMPLES
##################################

media_ee_c =NA
media_xx = NA
conta=  list(type = 'C0',ratio=0, media_xx = media_xx , media_ee_c =media_ee_c)
 
simulacion_ROC(Nrep=Nrep, nn = nn, grillat= grillat, modelo =1, conta=conta,
		   slope=slope, alfas = alfas, ee_sd =ee_sd,
                   criterio_seleccion = criterio_seleccion, funcionx0 = funcionx0,
		   rango_D=rango_D, rango_H=rango_H)
 
##############################
# CONTAMINATION C_{1,H}
##############################

media_ee_c = 8
media_xx = NA
media_vec =  media_ee_c*ee_sd
conta = list(type = 'C1H', ratio=0.1, xx_mu=NA, ee_c_mu = media_vec,ee_c_sd=0.5*ee_sd)


simulacion_ROC(Nrep=Nrep, nn = nn, grillat= grillat, modelo =1, conta=conta,
		   slope=slope, alfas = alfas, ee_sd =ee_sd,
                   criterio_seleccion = criterio_seleccion, funcionx0 = funcionx0,
		   rango_D=rango_D, rango_H=rango_H)


##############################
# CONTAMINATION C_{2,H}
##############################

media_ee_c = 12
media_xx = media_ee_c/2 
media_vec =  media_ee_c*ee_sd
conta = list(type = 'C2H', ratio=0.1,xx_col=2,xx_mu= media_vec/2 , xx_sd=0.5, ee_c_mu = media_vec,ee_c_sd=0.5*ee_sd)



simulacion_ROC(Nrep=Nrep, nn = nn, grillat= grillat, modelo =1, conta=conta,
		   slope=slope, alfas = alfas, ee_sd =ee_sd,
                   criterio_seleccion = criterio_seleccion, funcionx0 = funcionx0,
		   rango_D=rango_D, rango_H=rango_H)
