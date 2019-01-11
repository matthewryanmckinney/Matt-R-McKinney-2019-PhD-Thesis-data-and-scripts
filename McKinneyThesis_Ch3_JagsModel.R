model{
  for (j in 1:97){##nbrg
    for (t in 1:39){##nyears
      for (s in 1:9){##nsurvey types
      ###here y is ndetections in surveys
      
      y[j,t,s]~dbin(lo.psi[j,t,s],n[j,t,s]) ##n is nsurveys
      logit(lo.psi[j,t,s])<-psi[j,t,s]
      psi[j,t,s]~dnorm(mu.psi[j,t],tau.psi[t,s])
      
      }
    }
  }
  for (j in 1:97){
    for (t in 1:39){
      mu.psi[j,t]<-beta0[j] + beta1[j]*time[t] #+ beta2[j]*pow(time[t],2)
    }
  }
  for (t in 1:39){
    for (s in 1:9){
      tau.psi[t,s]~dunif(0,10)
    }
  }
  for (j in 1:97){###nbrgs
    beta0[j]~dnorm(mu.beta0[agstat[j]],tau.beta0[agstat[j]])
    beta1[j]~dnorm(mu.beta1[agstat[j]],tau.beta1[agstat[j]])
   # beta2[j]~dnorm(mu.beta2[agstat[j]],tau.beta2[agstat[j]])
   ##script for SSVS
    # beta2[j]~dnorm(0,b2.prec[j])
   # b2.prec[j]<-(gam2[j]*0.1)+((1-gam2[j])*99999999999999)
  #  gam2[j]~dbern(0.5)
    
   
    
  }
  for (b in 1:2){
    mu.beta0[b]~dnorm(mu.mu.beta0,tau.tau.beta0)
    mu.beta1[b]~dnorm(mu.mu.beta1,tau.tau.beta1)
    #mu.beta2[b]~dnorm(mu.mu.beta2,tau.tau.beta2)
    tau.beta0[b]~dgamma(0.01,0.01)
    tau.beta1[b]~dgamma(0.01,0.01)
    #tau.beta2[b]~dgamma(0.01,0.01)
  }
 
  mu.mu.beta0~dnorm(0,0.001)
  mu.mu.beta1~dnorm(0,0.001)
  #mu.mu.beta2~dnorm(0,0.01)
  tau.tau.beta0~dgamma(0.01,0.01)
  tau.tau.beta1~dgamma(0.01,0.01)
  #tau.tau.beta2~dgamma(0.01,0.01)
}##end model
  