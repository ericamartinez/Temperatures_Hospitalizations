
################################################################################
# "Temporal changes in the effects of ambient temperatures 
#               on hospital admissions in Spain"    
#   
# SECOND STAGE: multivariate meta-analysis  
#
#   ISGlobal  
#   October 2018
#
################################################################################


multivariate_ma <- function(dlist,coef1,vcov1,coef2,vcov2,provincies_n,varper,
                            prov_valides){
  
  ########################################################################################################
  # SECOND-STAGE ANALYSIS: MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF AND THEN COMPUTATION OF BLUP
  ########################################################################################################
  
  # CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
  avgtmean <- sapply(dlist,function(x) mean(x$tempmax_compl,na.rm=T))
  rangetmean <- sapply(dlist,function(x) diff(range(x$tempmax_compl,na.rm=T)))
  
  ################################################################################
  # META-ANALYSIS
  # We pool the estimated location-specific overall cumulative exposure-response associations 
  # (combina els efectes de cada provincia i obtenim una estimació global a nivell de país)
  
  # RUN THE MODELS FOR BEFORE AND AFTER PERIODS
  mv1 <- mvmeta(coef1~1,vcov1,provincies_n,control=list(showiter=T,maxiter=10000))
  summary(mv1)
  mv2 <- mvmeta(coef2~1,vcov2,provincies_n,control=list(showiter=T,maxiter=10000))
  summary(mv2)
  
  ################################################################################
  # OBTAIN BLUPS (best linear unbiased predictions)
  # For random effects models, predictions are the sum of the average results 
  # of the fixed part of the model plus the predicted random effects variances
  # We can draw the exposure-response curves in each city
  
  blup1 <- blup(mv1,vcov=T)
  blup2 <- blup(mv2,vcov=T)
  
  
  ################################################################################
  # RE-CENTERING
  
  # GENERATE THE MATRIX FOR STORING THE RESULTS
  minperccity <- mintempcity <- rep(NA,length(dlist))
  names(mintempcity) <- names(minperccity) <- provincies_n
  
  minperccity_1 <- minperccity_2 <- minperccity
  mintempcity_1 <- mintempcity_2 <- mintempcity
  
  predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
  tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                             predper/100,na.rm=T)))
  
  # DEFINE MINIMUM HOSPITALIZATIONS VALUES: EXCLUDE LOW AND VERY HOT TEMPERATURE
  ## BEFORE PERIOD ##
  for(i in seq(length(dlist))) {
    data <- dlist[[i]]
    data1 <- subset(data, data$yyyy<=2002)
    
    predvar <- quantile(data1$tempmax_compl,1:99/100,na.rm=T)
    # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
    
    # We define 3 knots at the 10th, 75th and 90th percentiles
    argvar <- list(x=predvar,fun="bs",
                   degree=2, knots=quantile(data1$tempmax_compl, varper/100, na.rm=T),
                   Bound=range(data1$tempmax_compl,na.rm=T))
    
    bvar <- do.call(onebasis,argvar)
    minperccity_1[i] <- (1:99)[which.min((bvar%*%blup1[[i]]$blup))]
    mintempcity_1[i] <- quantile(data1$tempmax_compl,minperccity_1[i]/100,na.rm=T)
  }
  
  # COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
  (minperccountry_1 <- median(minperccity_1))
  
  ## AFTER PERIOD ##
  for(i in seq(length(dlist))) {
    data <- dlist[[i]]
    data2 <- subset(data, data$yyyy>=2004)
    
    predvar <- quantile(data2$tempmax_compl,1:99/100,na.rm=T)
    # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
    
    # We define 3 knots at the 10th, 75th and 90th percentiles
    argvar <- list(x=predvar,fun="bs",
                   degree=2, knots=quantile(data2$tempmax_compl, varper/100, na.rm=T),
                   Bound=range(data2$tempmax_compl,na.rm=T))
    
    bvar <- do.call(onebasis,argvar)
    minperccity_2[i] <- (1:99)[which.min((bvar%*%blup2[[i]]$blup))]
    mintempcity_2[i] <- quantile(data2$tempmax_compl,minperccity_2[i]/100,na.rm=T)
  }
  
  # COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
  (minperccountry_2 <- median(minperccity_2))
  
  
  ################################################################################
  # PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
  
  ################################################################################
  
  prov <- c("Alava", "Albacete", "Alicante", "Almeria", "Avila", "Badajoz", "Illes Balears", "Barcelona", "Burgos",
            "Caceres", "Cadiz", "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", "Granada",
            "Guadalajara", "Guipuzkoa", "Huelva", "Huesca", "Jaen", "Leon", "Lleida", "La Rioja", "Lugo", "Madrid",
            "Malaga", "Murcia", "Navarra", "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
            "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
            "Valencia", "Valladolid", "Vizcaia", "Zamora", "Zaragoza")
  prov <- prov[prov_valides]
  
  datanew <- data.frame(avgtmean=mean(tapply(avgtmean,prov,mean)),
                        rangetmean=mean(tapply(rangetmean,prov,mean)))
  
  # PREDICT THE POOLED COEFFICIENTS FOR EACH MODEL
  mvpred1 <- predict(mv1,datanew,vcov=T,format="list")
  mvpred2 <- predict(mv2,datanew,vcov=T,format="list")
 
  
  # DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
  # ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
  predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
  set.seed(13041975)
  tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                             predper/100,na.rm=T)))
  
  # DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
  bvar <- onebasis(tmeancountry,fun="bs",degree=2, 
                   knots=tmeancountry[paste(varper,".0%",sep="")])
  
  aa=names(tmeancountry)[20:100]
  percmin1=aa[which.min((bvar%*%mvpred1$fit)[20:100])]
  tempmin1=tmeancountry[percmin1]
  percmin2=aa[which.min((bvar%*%mvpred2$fit)[20:100])]
  tempmin2=tmeancountry[percmin2]
  
  #cenindcountry1 <- which.min(bvar%*%mvpred1$fit)
  #cenindcountry2 <- which.min(bvar%*%mvpred2$fit)
  
  # DEFINE CENTERING PERCENTILE FOR COUNTRY
  #cenpercountry1 <- pmin(pmax(predper[cenindcountry1],10),90)
  #cenpercountry2 <- pmin(pmax(predper[cenindcountry2],10),90)
  
  
  ################################################################################
  # PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
  
  # OBTAIN THE CENTERED PREDICTIONS
  #bvar <- onebasis(tmeancountry,fun="bs",degree=2,
  #                 knots=tmeancountry[paste(varper,".0%",sep="")])
 
  #cen1 <- tmeancountry[paste(cenpercountry1,".0%",sep="")]
  #cen2 <- tmeancountry[paste(cenpercountry2,".0%",sep="")]
  
  cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,model.link="log",
                   at=tmeancountry,cen=tempmin1)
  cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,model.link="log",
                   at=tmeancountry,cen=tempmin2)
  
  
  return(list(cp1,cp2,blup1,blup2,tmeancountry,percmin1,percmin2))
  
  #
}


