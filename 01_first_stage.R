
################################################################################
# "Temporal changes in the effects of ambient temperatures 
#               on hospital admissions in Spain"    
#
#   FIRST STAGE: time series models for each province
#   
#   ISGlobal  
#   October 2018
#
################################################################################

tsmodels_provinces <- function (dlist,provinces,provincies_n,varper,vardegree,lag,arglag,
                         name_var,dftrend,dfseas,j,variables_fallen){
  
  ################################################################################
  # FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH PROVINCE, REDUCE AND SAVE
  ################################################################################
  
  ################################################################################
  # CREATE THE OBJECTS TO STORE THE RESULTS
  
  # COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
  coef <- coef1 <- coef2 <- coefint <- matrix(NA,length(provinces), length(varper)+vardegree,
                                              dimnames= list(provincies_n))
  
  # The matrix has 5 columns
  
  vcov <- vcov1 <- vcov2 <- vcovint <- vector("list",length(provinces))
  names(vcov) <- names(vcov1) <- names(vcov2) <- names(vcovint) <- provincies_n
  
  
  ################################################################################
  # RUN THE LOOP
  
  # LOOP
  time <- proc.time()[3]
  for(i in seq(length(dlist))) {
    try({
    # PRINT
    cat(i,"")
    
    # EXTRACT THE DATA
    data <- dlist[[i]]
      
    # BEFORE PERIOD (1997-2002)
    data1 <- subset(data, data$yyyy<=2002)
    eval(parse(text = paste0("outcome1 <- data1$", name_var)))
    data1$t=1:dim(data1)[1]
    
    argvar1 <- list(fun="bs",degree=2,knots=quantile(data1$tempmax_compl,varper/100,na.rm=T))
    
    cb1 <- crossbasis(data1$tempmax_compl, lag=lag, 
                      argvar=argvar1,arglag=arglag)
    
    model1 <- glm(outcome1 ~ cb1 + dow + hday + phday + total_influenza_h +
                    ns (doy,df=dfseas): factor(yyyy) + 
                    ns(t, df=round(length(unique(data1$yyyy))/dftrend/10)), 
                  data1, family=quasipoisson, na.action="na.exclude",
                  control = list(maxit = 5000))
    
    # AFTER PERIOD (2004-2013)
    data2 <- subset(data, data$yyyy>=2004)
    eval(parse(text = paste0("outcome2 <- data2$", name_var)))
    data2$t=1:dim(data2)[1]
    
    argvar2 <- list(fun="bs",degree=2,knots=quantile(data2$tempmax_compl,varper/100,na.rm=T))
    
    cb2 <- crossbasis(data2$tempmax_compl, lag=lag, 
                      argvar=argvar2, arglag=arglag)
    
    model2 <- glm(outcome2 ~ cb2 + dow + hday + phday + total_influenza_h +
                    ns (doy,df=dfseas): factor(yyyy) + 
                    ns(t, df=round(length(unique(data2$yyyy))/dftrend/10)), 
                  data2, family=quasipoisson, na.action="na.exclude",
                  control = list(maxit = 5000))
    
    
    # REDUCTION TO OVERALL CUMULATIVE
    # Sum the effects of all lags in order to eliminate one dimension of the association. 
    #Sum (acumulate) the risk during the lag period
    
    red1 <- crossreduce(cb1,model1)
    coef1[i,] <- coef(red1)
    vcov1[[i]] <- vcov(red1)
    
    red2 <- crossreduce(cb2,model2)
    coef2[i,] <- coef(red2)
    vcov2[[i]] <- vcov(red2)
    }) 
  }
  proc.time()[3]-time
  
  prov_valides <- !is.na(coef1[,1])
  dlist_valides <- dlist[prov_valides]
  provincies_n_valides <- provincies_n[prov_valides]
  coef1_valides <- coef1[prov_valides,]
  coef2_valides <- coef2[prov_valides,]
  vcov1_valides <- vcov1[prov_valides]
  vcov2_valides <- vcov2[prov_valides]
  provinces_valides <- provinces[prov_valides]
  variables_fallen[j,1]  <- name_var
  variables_fallen[j,2:51] <- prov_valides
  
  return(list(coef1_valides,vcov1_valides,coef2_valides,
              vcov2_valides,prov_valides,
              dlist_valides,provincies_n_valides,provinces_valides,variables_fallen))
  #
}

