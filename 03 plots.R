
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
# PLOTS
################################################################################


# POOLED: BEFORE AND AFTER PERIODS

period <- c("Period 1","Period 2")
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
indlab <- predper%in%c(0,1,10,50,90,99,100)

plot(cp1,ylab="Percent change (%)",xlab="Temperature percentile",axes=F,
     ylim=c(0.6,1.8),lwd=2,col="chocolate1",ci.arg=list(density=20,col="chocolate1"))
lines(cp2,ci="area",lwd=2,col="olivedrab4",ci.arg=list(density=20,angle=-45,col="olivedrab4"))
axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9, at=c(0.6,0.7,0.8,0.9,1.00,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8), 
     labels=c(-40,-30,-20,-10,0,10,20,30,40,50,60,70,80))
abline(v=tmeancountry[c("1.0%","99.0%")],
       lty=c(2,2))
legend("top",period,col=c("chocolate1","olivedrab4"),lwd=2,bg="white",cex=0.8,ncol=2)



####################
##
## 1. MAP NUMBER OF DAYS THAT THE PLAN WAS ACTIVATED####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm

library("sp")

esp1 <- readRDS("ESP_adm2.rds")

esp1@data$nump1 <- c(70,11,15,0,25,5,88,15,69,0,11,11,0,39,131,36,0,71,0,0,5,0,0,
                     24,0,5,0,21,10,5,10,0,0,5,0,10,5,20,25,20,5,0,0,0,10,10,10,10,5,5,0,0,0,0,0,20)
esp1@data$nump2 <- c(158,118,68,5,86,37,169,119,197,0,90,90,0,121,240,81,17,120,0,5,5,5,15,94,0,5,5,223,73,40,
                     25,0,0,65,10,15,56,35,86,96,5,0,0,11,31,74,110,110,6,6,5,5,0,5,0,145)

# Period 1
levels <- cut(esp1@data$nump1,c(-50,10,50,300),
              labels=c('<10', '[10-50]', '>50'))
funcol <- colorRampPalette(c("chartreuse3", "darkorange","red2"))
col <- funcol(length(levels(levels)))[levels]

plot(esp1, col=col)
text(coordinates(esp1), labels = esp1@data$nump1, cex=0.3)
legend(-12,39,paste0(c('<10', '[10-50]', '>50')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Number of days")	


# Period 2
levels <- cut(esp1@data$nump2,c(-50,10,50,300),
              labels=c('<10', '[10-50]', '>50'))
funcol <- colorRampPalette(c("chartreuse3", "darkorange","red2"))
col <- funcol(length(levels(levels)))[levels]

plot(esp1, col=col)
text(coordinates(esp1), labels = esp1@data$nump2, cex=0.3)
legend(-12,39,paste0(c('<10', '[10-50]', '>50')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Number of days")	
