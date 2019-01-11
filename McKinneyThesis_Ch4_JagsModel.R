model{
  for (s in 1:2357){##surveys
    for (i in 1:Nspecies){##species; spp 1 is Myna, Spp2 is Miner, #3 is magpie lark, #4 is rainbow lorikeet
      y[s,i]~dbin(p[s,i],z[s,i])
      z[s,i]~dpois(eff.lambda[site[s],i,year[s]])
      logit(p[s,i])<-alpha1[type[s],i] +
        #alpha2[time[s],i] +
        alpha4[i]*date[s] + alpha5[i]*date[s]*date[s]
    }
  }
  for (j in 1:30){##sites
    for (i in 1:Nspecies){##species
      for (t in 1:24){##years
        eff.lambda[j,i,t]<-lambda[j,i,t]*wz[j,i,t]
        wz[j,i,t]~dbern(psi[i,t])
        log(lambda[j,i,t])<-n[j,i,t]
        
        ######## Incorporating demographic stochasticity ################# 
        n[j,i,t]~dnorm(w[j,i,t], Dem.prec[j,i,t])
        Dem.var[j,i,t]<-1/Dem.prec[j,i,t]
        
      }
    }
  }
  for (i in 1:Nspecies){
    for (t in 1:24){
      psi[i,t]~dunif(0,1)
    }
  }
  for (j in 1:30){
    for (t in 1:24){
      w[j,1:Nspecies,t]~dmnorm(m[j,1:Nspecies,t], Tau[,])
    }
  }
  for (j in 1:30){##sites
    for (i in 1:Nspecies){##species
      for (t in 2:24){##years
        Dem.prec[j,i,t]<-tau.d[i]*exp(n[j,i,t-1]) ###tau.d is delta in publication
        m[j,i,t]<- n[j,i,t-1] + r[i]*(1-(inprod(alpha[i,], n[j,,t-1])/k[i]))+ 
          beta1[i]*dCoast[j] +
          beta2[i]*lat[j] #+
        #beta3[i]*imp[j] 
      }
    }
  }
  for (j in 1:30){##sites
    for (i in 1:Nspecies){##species
      m[j,i,1]~dnorm(0,0.001)
      Dem.prec[j,i,1]~dunif(0,15)
    }
  }
  ############## PRIOR ELICITATION ####################################### 
  for (i in 1:Nspecies){
    ### Proportions of variance explained ################################# 
    cenv[i]<-Sigma[i,i]/totvar[i]
    cintra[i]<-intra[i]/totvar[i]
    cinter[i]<-inter[i]/totvar[i]
    ### Total variance for species i #######################################
    totvar[i]<-pow((r[i]/k[i]),2)*inprod(pow(alpha[i,],2),svar[])+ct[i] + Sigma[i,i]
    ## svar[i] is the MCMC-based (stationary) variance of n[,i] from a preliminary analysis
    ### Prop of env. Var. explained by the covariates ####################
    p.m.env1[i]<-(beta1[i]*beta1[i])/env[i]
    p.m.env2[i]<-(beta2[i]*beta2[i])/env[i]
    #p.m.env3[i]<-(beta3[i]*beta3[i])/env[i]
    
    p.over.env[i]<-ct[i]/env[i] # jointly
    
    #### Env. stochast. and intra/inter-spec. interactions ############### 
    env[i]<-ct[i]+Sigma[i,i]
    ######## Covariance explained by the environmental surrogate covariates ############
    ct[i]<-beta1[i]*beta1[i] + beta2[i]*beta2[i] #+ beta3[i]*beta3[i]
    
    inter[i]<-totinteract[i]-intra[i]
    intra[i]<-pow((r[i]/k[i]),2)*svar[i]
    totinteract[i]<-pow((r[i]/k[i]),2)*inprod(pow(alpha[i,],2),svar[])
  }
  ### Prior for the residual environmental covariance matrix ########## 
  Sigma[1:Nspecies,1:Nspecies]<-inverse(Tau[,])
  Tau[1:Nspecies,1:Nspecies] ~ dwish(Rho[ , ], Nspecies)##Rho is an identity matrix specified in the data
  ## priors for interaction coefficients ################################ 
  for(i in 1:Nspecies){ 
    alpha[i,i]<-1
    gam[i,i]<-1
  }
  for(i in 1:(Nspecies-1)){
    ### gam is the inclusion indicator for SSVS ################################# 
    for(j in (i+1):Nspecies){  
      alpha[i,j]~dnorm(0, prec[i,j])
      alpha[j,i]~dnorm(0, prec[j,i])
      prec[i,j]<-1/vari[i,j] 
      prec[j,i]<-1/vari[j,i] 
      vari[i,j]<-(1-gam[i,j])*0.00000000000001+gam[i,j]*10 
      vari[j,i]<-(1-gam[j,i])*0.00000000000001+gam[j,i]*10
      gam[i,j]~dbern(0.2) 
      gam[j,i]~dbern(0.2)
    } 
  }
  #alpha[1,2]~dnorm(0, prec[1,2])
  #alpha[2,1]~dnorm(0, prec[2,1])
  #prec[1,2]<-1/vari[1,2] 
  #prec[2,1]<-1/vari[2,1] 
  #vari[1,2]<-((1-gam[1,2])*0.00000000000001)+gam[1,2]*10 
  #vari[2,1]<-((1-gam[2,1])*0.00000000000001)+gam[2,1]*10
  #gam[1,2]~dbern(0.2) 
  #gam[2,1]~dbern(0.2)
  k[1]~dunif(3,5.5)
  k[2]~dunif(3,4.5)
  k[3]~dunif(3,6.8)
  k[4]~dunif(3,4.5)
  ### Priors for other variables ########################################## 
  
  for(i in 1:Nspecies){
    r[i]~dunif(0, 2)
    
    tau.d[i]<-1/vd[i]
    vd[i]<-sdem[i]*sdem[i]
    sdem[i]~dunif(0.01, 10)
    
    alpha0[i]~dnorm(0,0.01)
    beta1[i]~dnorm(0,prec.b1[i])
    beta2[i]~dnorm(0,prec.b2[i])
    
    prec.b1[i]<-1/vari.b1[i]
    prec.b2[i]<-1/vari.b2[i]
    
    vari.b1[i]<-((1-gam.b1[i])*0.00000000000001)+gam.b1[i]*10 
    vari.b2[i]<-((1-gam.b2[i])*0.00000000000001)+gam.b2[i]*10
    gam.b1[i]~dbern(0.2) 
    gam.b2[i]~dbern(0.2)
    
    
    #beta3[i]~dnorm(0,0.01)
    alpha4[i]~dnorm(0,0.01)
    alpha5[i]~dnorm(0,0.01)
  }
  for (f in 1:ntype){##6 types
    for (i in 1:Nspecies){
      alpha1[f,i]~dnorm(0,0.01)
    }
  }
  for (f in 1:ntime){##8 times
    for (i in 1:Nspecies){
      alpha2[f,i]~dnorm(0,0.01)
    }
  }
}##end model