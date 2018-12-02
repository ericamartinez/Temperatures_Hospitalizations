
################################################################################
# "Temporal changes in the effects of ambient temperatures 
#               on hospital admissions in Spain"    
#
#   ISGlobal  
#   October 2018
#   
#
################################################################################


################################################################################
# PREPARE THE DATA
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel);library(mgcv)
library(foreach); library(doSNOW); library(Cairo); library(parallel)

load("Hospit_temp.Rdata")

# Day of the year (1:365)
hospit_temp$date <- paste(hospit_temp$dd, hospit_temp$mm, hospit_temp$yyyy, sep="/")
hospit_temp$date <- strptime(hospit_temp$date, '%d/%m/%Y')
hospit_temp$doy <- as.numeric(strftime(hospit_temp$date, format = "%j"))

# We exclude 2003
hospit_temp <- subset(hospit_temp, hospit_temp$yyyy!=2003)

hospit_temp$dow <- as.factor(hospit_temp$dow)
hospit_temp <- hospit_temp[with(hospit_temp, order(provcode, yyyy, mm, dd)), ]

###Pick holidays (summer and Christmas)
hospit_temp$phday <- 0

# First three weeks of August
hospit_temp[hospit_temp$dd < 22 & hospit_temp$mm == 8, 'phday'] <- 1
# From the 23rd of December until the 6th of January
hospit_temp[(hospit_temp$dd > 22 & hospit_temp$mm == 12) | 
              (hospit_temp$dd < 7 & hospit_temp$mm == 1), 'phday'] <- 1

# Binary variables for before and after periods
hospit_temp$int1 <- 0
hospit_temp[hospit_temp$yyyy >= 2004, 'int1'] <- 1

hospit_temp$int2 <- 0
hospit_temp[hospit_temp$yyyy < 2004, 'int2'] <- 1


provincies_n_total <- list("Alava", "Albacete", "Alicante", "Almeria", "Avila", 
                           "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
                           "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
                           "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
                           "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
                           "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
                           "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
                           "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
                           "Zamora", "Zaragoza")


# ARRANGE THE DATA AS A LIST OF DATA SETS
provinces_total <- as.character(unique(hospit_temp$provcode)) # My provinces
dlist_total <- lapply(provinces_total,function(x) hospit_temp[hospit_temp$provcode==x,]) 
# Create a list with 50 provinces 
#(agafa el data frame i el converteix a llista de tants elements com provincies)
names(dlist_total) <- provincies_n_total

#### PARAMETERS FOR THE MAIN MODEL
# 3 internal knots placed at the 10th, 50th and 90th percentiles of 
# location-specific temperature distribution
varper <- c(10,50,90)
vardegree <- 2


# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3 #Number of knots for lag model
arglag<- list(knots=logknots(lag,lagnk))

# DEEGRES OF FREEDOM
dfseas <- 8  #Seasonality
dftrend <- 1 #Long-term trend



# LOAD THE FUNCTION FOR COMPUTING THE ATTRIBUTABLE RISK MEASURES
source("01_first_stage.R")
source("02_second_stage.R")


variables_fallen_total <- matrix(NA,nrow=195,ncol=51)
colnames(variables_fallen_total) <- c("Variable","Alava", "Albacete", "Alicante", "Almeria", "Avila", 
                                      "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
                                      "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
                                      "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
                                      "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
                                      "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
                                      "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
                                      "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
                                      "Zamora", "Zaragoza")

variables_article <- c("total_cvd_h","a1664cvd_h","a6574cvd_h","a7584cvd_h","a85cvd_h",   #cardiovascular
                     "total_cere_h","a1664cere_h","a6574cere_h","a7584cere_h","a85cere_h",   #cerebrovascular
                     "total_resp_h","a1664resp_h","a6574resp_h","a7584resp_h","a85resp_h")   #respiratory   

variables_final_total <- matrix(NA,nrow=195,ncol=51)
colnames(variables_final_total) <- c("Variable","Alava", "Albacete", "Alicante", "Almeria", "Avila", 
                                      "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
                                      "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
                                      "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
                                      "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
                                      "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
                                      "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
                                      "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
                                      "Zamora", "Zaragoza")



#####################################################################################
##
## RESULTS VARIABLES ####
##
#####################################################################################

results <- matrix (NA, nrow=length(variables_article), ncol=15)
colnames(results) <- c("variable","RR-1st P1", "IClow - 1st P1", "ICup - 1st P1", "RR-99th P1", 
                       "IClow - 99th P1", "ICup - 99th P1", "RR-1st P2", "IClow - 1st P2", 
                       "ICup - 1st P2", "RR-99th P2", "IClow - 99th P2", "ICup - 99th P2", "MMP1","MMP2")

results_moderate <- matrix (NA, nrow=length(variables_article), ncol=14)
colnames(results_moderate) <- c("variable","RR-10th P1", "IClow - 10th P1", "ICup - 10th P1", "RR-90th P1", 
                       "IClow - 90th P1", "ICup - 90th P1", "RR-10th P2", "IClow - 10th P2", 
                       "ICup - 10th P2", "RR-90th P2", "IClow - 90th P2", "ICup - 909th P2","Prov_valides")

for (j in 1:10){ 
  
  name_var <- variables_article[j]
  print(variables_article[j])
  
  ### FIRST STAGE ###
  ts <- tsmodels_provinces(dlist_total,provinces_total,provincies_n_total,varper,vardegree,lag,arglag,
                           name_var,dftrend,dfseas,j,variables_final_total)
  coef1 <- ts[[1]]; vcov1 <- ts[[2]]; coef2 <- ts[[3]]
  vcov2 <- ts[[4]]; prov_valides <- ts[[5]]; dlist <- ts[[6]]; 
  provincies_n <- ts[[7]]; provinces <- ts[[8]]; variables_final <- ts[[9]]
  
  print(sum(prov_valides))
  
  ### SECOND STAGE ###
  mv <- multivariate_ma(dlist,coef1,vcov1,coef2,vcov2,
                        provincies_n,varper,prov_valides)
  cp1 <- mv[[1]]; cp2 <- mv[[2]]; 
  blup1 <- mv[[3]]; blup2 <- mv[[4]];
  tmeancountry <- mv[[5]]; 
  cenpercountry1 <- mv[[6]]; cenpercountry2 <- mv[[7]]
  
  ### MAIN RESULTS ###
  predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
  results[j,1] <- name_var
  # RR AT 1st AND 99TH VS MMP (WITH 95%CI) - EXTREME TEMPERATURES
  predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
  #### BEFORE PERIOD ####
  results[j,2] <- cp1$allRRfit[predper==1];results[j,3] <- cp1$allRRlow[predper==1];results[j,4] <- cp1$allRRhigh[predper==1]
  results[j,5] <- cp1$allRRfit[predper==99];results[j,6] <- cp1$allRRlow[predper==99];results[j,7] <- cp1$allRRhigh[predper==99]
  #### AFTER PERIOD ####
  results[j,8] <- cp2$allRRfit[predper==1];results[j,9] <- cp2$allRRlow[predper==1];results[j,10] <- cp2$allRRhigh[predper==1]
  results[j,11] <- cp2$allRRfit[predper==99];results[j,12] <- cp2$allRRlow[predper==99];results[j,13] <- cp2$allRRhigh[predper==99]
  # Minimum hospitalization percentile
  results[j,14] <- cenpercountry1
  results[j,15] <- cenpercountry2
  
  # RR AT 10th AND 90TH VS MMP (WITH 95%CI) - MODERATE TEMPERATURES
  results_moderate[j,1] <- name_var
  #### BEFORE PERIOD ####
  results_moderate[j,2] <- cp1$allRRfit[predper==10];results_moderate[j,3] <- cp1$allRRlow[predper==10];results_moderate[j,4] <- cp1$allRRhigh[predper==10]
  results_moderate[j,5] <- cp1$allRRfit[predper==90];results_moderate[j,6] <- cp1$allRRlow[predper==90];results_moderate[j,7] <- cp1$allRRhigh[predper==90]
  #### AFTER PERIOD ####
  results_moderate[j,8] <- cp2$allRRfit[predper==10];results_moderate[j,9] <- cp2$allRRlow[predper==10];results_moderate[j,10] <- cp2$allRRhigh[predper==10]
  results_moderate[j,11] <- cp2$allRRfit[predper==90];results_moderate[j,12] <- cp2$allRRlow[predper==90];results_moderate[j,13] <- cp2$allRRhigh[predper==90]
  
  results_moderate[j,14] <- sum(prov_valides)
  
}



