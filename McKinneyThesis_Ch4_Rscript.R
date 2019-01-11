
##McKinney 2019 Thesis Chapter 3 script for running analysis
##data and model can be found at 
#y4 is matrix of individuals observed for each survey (rows) for each species (columns)
#site is indexed site ID
#year is indexed year ID
#type is survey type
#time is time amount for that survey
#date is julian date for survey, standardized to mean =0 and sd=1
#dCoast is distance to coast for each site, standardized
#lat is latitude of each site, standardized
#ntime is the number of different time variables
#nSpecies4 = 4, the number of species in this analysis
#ntype is the number of data types
#rho4 is the 4x4 identity matrix
#svar4 is the stationary variance from a preliminary analysis

library(snow)
library(rjags)
library(runjags)
library(DHARMa)
library(R2jags)


forjags<-list("y"=y4,"site"=site,"year"=year,"type"=type,"time"=time,"date"=date,
              "dCoast"=dCoast,"lat"=lat,"ntime"=ntime,"Nspecies"=Nspecies4,"ntype"=ntype,"Rho"=rho4,"svar"=svar4)

b.ins<-array(rnorm((4*3*3),0,0.01),dim=c(4,3,3))
a1.ins<-array(rnorm((4*6*3),-1,0.01),dim=c(6,4,3))
a2.ins<-array(rnorm((4*8*3),-1,0.01),dim=c(8,4,3))
a4.ins<-array(rnorm((4*3),0,0.01),dim=c(4,3))
a5.ins<-array(rnorm((4*3),0,0.01),dim=c(4,3))
wz.ins<-array(rep(1,(30*4*24)),dim=c(30,4,24))
z.ins<-z.ins4
gam.ins<-matrix(data = rep(1,4),nrow = 2,ncol = 2)
d.ins<-matrix(data=rep(1,2357*2),nrow=2357,ncol=2)



j.inits<-list(list(beta1=b.ins[,1,1],beta2=b.ins[,2,1],beta3=b.ins[,3,1],
                   alpha1=a1.ins[,,1],alpha2=a2.ins[,,1],alpha4=a4.ins[,1],alpha5=a5.ins[,1],
                   wz=wz.ins,z=((y4)+1),d=d.ins),
              list(beta1=b.ins[,1,2],beta2=b.ins[,2,2],beta3=b.ins[,3,2],
                   alpha1=a1.ins[,,2],alpha2=a2.ins[,,2],alpha4=a4.ins[,2],alpha5=a5.ins[,2],
                   wz=wz.ins,z=((y4)+1),d=d.ins),
              list(beta1=b.ins[,1,3],beta2=b.ins[,2,3],beta3=b.ins[,3,3],
                   alpha1=a1.ins[,,3],alpha2=a2.ins[,,3],alpha4=a4.ins[,3],alpha5=a5.ins[,3],
                   wz=wz.ins,z=((y4)+1),d=d.ins))

parms<-c("cenv","cintra","cinter","totvar","p.m.env1","p.m.env2","p.over.env","env",
         "ct","inter","intra","totinteract","alpha","gam","r","k","beta1","gam.b1","beta2","gam.b2")
memory.limit(size=10000000)
cjagsout<-run.jags(model="C:/Users/uqmmckin/Documents/Code/MynaMiner05122018_4spp_z.R",
                   monitor=parms,data=forjags,n.chains=3,burnin=500000,adapt=1000,inits=j.inits,
                   sample=5000,thin=10,jags.refresh=60,method="rjags",modules = "glm on",summarise = TRUE)


results<-summary(cjagsout)
results<-data.frame(results)



