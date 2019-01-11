##McKinney 2019 Thesis Chapter 3 script for running analysis
##data can be found at 
#ndet.XXXX1980 is the number of detections for: Star = starlings, Hspar=House Sparrows, Tspar=Tree Sparrows, Myna = Mynas
#time is the time in years elapsed since 1980, standardized to mean 0 and SD=1
#nsurv1980 is the number of surveys conducted
#agstat1980 is the agricultural status (ag=2, non-ag=1) for each bioregion
#
library(snow)
library(rjags)
library(runjags)
library(DHARMa)
library(R2jags)



forjags<-list("y"=ndet.Star1980,"time"=time.y1980,"n"=nsurv1980,"agstat"=agstat1980)
b.ins<-matrix(rnorm((97*3*3),0,0.01),nrow=97,ncol=9)

#z.ins<-as.numeric(c(rep(1,599777)))
j.inits<-list(list(beta0=b.ins[,1],beta1=b.ins[,2],beta2=b.ins[,3]),
              list(beta0=b.ins[,4],beta1=b.ins[,5],beta2=b.ins[,6]),
              list(beta0=b.ins[,7],beta1=b.ins[,8],beta2=b.ins[,9]))

parms<-c("beta0","beta1","mu.beta0","mu.beta1")
memory.limit(size=10000000)
cjagsout<-run.jags(model="C:/Users/uqmmckin/Documents/Code/CavNest_trend_simp_red.R",
                   monitor=parms,data=forjags,n.chains=3,burnin=100000,adapt=1000,inits=j.inits,
                   sample=5000,thin=10,jags.refresh=60,method="rjags",modules = "glm on",summarise = TRUE)


results<-summary(cjagsout)
