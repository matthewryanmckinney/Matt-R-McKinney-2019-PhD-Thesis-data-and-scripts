###McKinney Thesis 2019; Chapter 2 analysis
library(hier.part)
library(boot)
library(Hmisc)
library(MASS)
library(nlme)
#Data can be found at
#NDVI_mean is the mean NDVI for the Australia-wide analysis
#temp is temperature for the Australia-wide analysis
#hfp is human footprint for Australia-wide analysis
#nativerich is the native bird species richness for the Australia-wide analysis
#island is categorical variable indicating island or not for the Australia-wide analysis
#intros is the number of introduction events for the Australia-wide analysis
#HabDiv is the landcover type diversity for the Australia-wide analysis
#state is the categorical variable of Australian state for the Australia-wide analysis
#alienrich is the alien bird species richness for tthe Australia-wide analysis
#InGen is a data.frame that holds all of the above variables
#all of the above objects have an associated object with a "T" suffix for the Tasmania-only analysis
######################################################################
#Hierarchical partitioning
######################################################################
#General;Entire Continent
multix<-data.frame(NDVI_Mean,temp,hfp,nativerich,island,intros,HabDiv,state)

names(multix)
y<-(alienrich)
hp<-hier.part(y,multix,family="poisson",gof="logLik")

summary(hp)
hp$I.perc
hpgfs<-all.regs(y,multix,family="poisson",gof="logLik")
partition(hpgfs, pcan = 8, var.names = names(multix))

#General;Tazzy Only;
multixT<-data.frame(NDVI_Meant,tempt,hfpt,nativericht,islandt,introst,HabDivt)

names(multixT)
yt<-(alienricht)
hpT<-hier.part(yt,multixT,family="poisson",gof="logLik")

hpT$I.perc

####################################################################################################
#Randomization tests!!!This takes a LONG TIME!!! 
#Bootstrap confidence intervals!!!This takes a LONG TIME!!! 
#Does not work yet
#####################################################################################################
randhp<-rand.hp(y,multix,family="poisson",gof="logLik",num.rep=1000)
randhp$Iprobs

randhpTazzy<-rand.hp(yt,multixT,family="poisson",gof="logLik",num.rep=1000)
randhpTazzy$Iprobs


###GLS to test for spatial covariance structure

########################################################################
##Models for general analysis, Australia-wide
############################################################################
##Non-spatial model
glsIG<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv+State,data=InGen)
AIC(glsIG)
##spatial model--exponential covariance
glsIGexp<-update(glsIG,correlation=corExp(1,form=~POINT_X+POINT_Y))
AIC(glsIGexp)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
##spatial model--Gaussian covariance
glsIGgau<-update(glsIG,correlation=corGaus(1,form=~POINT_X+POINT_Y))
AIC(glsIGgau)
##Spatial model--spherical covariance
glsIGsph<-update(glsIG,correlation=corSpher(1,form=~POINT_X+POINT_Y))
AIC(glsIGsph)
##spatial model--linear covariance
glsIGlin<-update(glsIG,correlation=corLin(1,form=~POINT_X+POINT_Y))
AIC(glsIGlin)
#plot(Variogram(glsIG,form=~POINT_X+POINT_Y))
########################################################################
##Models for general analysis,Tazzy-only
############################################################################
##Non-spatial model
glsIGT<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,data=InGenT)
AIC(glsIGT)
##spatial model--exponential covariance
glsIGTexp<-update(glsIGT,correlation=corExp(1,form=~POINT_X+POINT_Y))
AIC(glsIGTexp)
##spatial model--Gaussian covariance
glsIGTgau<-update(glsIGT,correlation=corGaus(1,form=~POINT_X+POINT_Y))
AIC(glsIGTgau)
##Spatial model--spherical covariance
glsIGTsph<-update(glsIGT,correlation=corSpher(1,form=~POINT_X+POINT_Y))
AIC(glsIGTsph)
##spatial model--linear covariance
glsIGTlin<-update(glsIGT,correlation=corLin(1,form=~POINT_X+POINT_Y))
AIC(glsIGTlin)

########################################################################
##AIC/BIC for models
#########################################################################
#General,Australia wide
AIC(glsIG,glsIGexp,glsIGgau,glsIGsph)
BIC(glsIG,glsIGexp,glsIGgau,glsIGsph,glsIGlin)
#General, TAZZY Only
AIC(glsIGT,glsIGTexp,glsIGTgau,glsIGTsph)
BIC(glsIGT,glsIGTexp,glsIGTgau,glsIGTsph)

############################################################################
############################################################################
##Correlations--model residuals vs. excluded variables
############################################################################
############################################################################
##General Analysis, Australia-wide
############################################################################
##spatial model--exponential covariance; NDVI_Mean EXCLUDED
glsIGexp1<-gls(AlienSppCounts_post2000_50km~MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr1<-rcorr(InGen$NDVI_Mean,glsIGexp1$resid,type="spearman")
pr1
##spatial model--exponential covariance; MinTemp EXCLUDED
glsIGexp2<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr2<-rcorr(InGen$MinTemp,glsIGexp2$resid,type="spearman")
pr2
##spatial model--exponential covariance; hfp_clip EXCLUDED
glsIGexp3<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr3<-rcorr(InGen$hfp_clip50,glsIGexp3$resid,type="spearman")
pr3
##spatial model--exponential covariance; NativeSppCounts2014_50k EXCLUDED
glsIGexp4<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+hfp_clip50+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr4<-rcorr(InGen$NativeSppCounts2014_50k,glsIGexp4$resid,type="spearman")
pr4
##spatial model--exponential covariance; Island_Main EXCLUDED
glsIGexp5<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+hfp_clip50+NativeSppCounts2014_50k+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr5<-rcorr(InGen$Island_Main,glsIGexp5$resid,type="spearman")
pr5
##spatial model--exponential covariance; SppIntroCountsAliens50k EXCLUDED
glsIGexp6<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr6<-rcorr(InGen$SppIntroCountsAliens50k,glsIGexp6$resid,type="spearman")
pr6
##spatial model--exponential covariance; HabDiv EXCLUDED
glsIGexp7<-gls(AlienSppCounts_post2000_50km~NDVI_Mean+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGen)

pr7<-rcorr(InGen$HabDiv,glsIGexp7$resid,type="spearman")
pr7
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
############################################################################
##General Analysis, Tazzy-only
############################################################################
##Spatial model, exponential covariance; NDVI_Mean EXCLUDED
glsIGT1<-gls(AlienSppCounts_post2000_50km~MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr8<-rcorr(InGenT$NDVI_Mean50,glsIGT1$resid,type="spearman")
pr8
##Spatial model, exponential covariance; MinTemp EXCLUDED
glsIGT2<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+hfp_clip50+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr9<-rcorr(InGenT$MinTemp,glsIGT2$resid,type="spearman")
pr9
##Spatial model, exponential covariance; hfp_clip EXCLUDED
glsIGT3<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+NativeSppCounts2014_50k+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr10<-rcorr(InGenT$hfp_clip50,glsIGT3$resid,type="spearman")
pr10
##Spatial model, exponential covariance; NativeSppCounts2014_50k EXCLUDED
glsIGT4<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+hfp_clip50+Island_Main+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr11<-rcorr(InGenT$NativeSppCounts2014_50k,glsIGT4$resid,type="spearman")
pr11
##Spatial model, exponential covariance; Island_Main EXCLUDED
glsIGT5<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+hfp_clip50+NativeSppCounts2014_50k+SppIntroCountsAliens50k+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr12<-rcorr(InGenT$Island_Main,glsIGT5$resid,type="spearman")
pr12
##Spatial model, exponential covariance; SppIntroCountsAliens50k EXCLUDED
glsIGT6<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main+HabDiv,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr13<-rcorr(InGenT$SppIntroCountsAliens50k,glsIGT6$resid,type="spearman")
pr13

##Spatial model, exponential covariance; HabDiv EXCLUDED
glsIGT7<-gls(AlienSppCounts_post2000_50km~NDVI_Mean50+MinTemp+hfp_clip50+NativeSppCounts2014_50k+Island_Main,correlation=corExp(1,form=~POINT_X+POINT_Y),data=InGenT)
#save.image("~/Rworkspaces/HPGLS50kEuroAus_ProcB.RData")
pr14<-rcorr(InGenT$HabDiv,glsIGT7$resid,type="spearman")
pr14
