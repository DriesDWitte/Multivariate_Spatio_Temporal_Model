################################################################################
# This script contains the code to fit the multivariate spatio-temporal model  #
# using the MCAR prior for the spatial random effects and the bivariate normal #
# prior for temporal random effects.                                           #                   #
################################################################################

#load packages
library(sf)
library(rgdal)
library(ggplot2)
library(here)
library(readxl)
library(dplyr)
library(lubridate)
library(tidyr)
library(CARBayes)
library(CARBayesST)
library(GGally)
library(spdep)
library (readr)
library(coda)
library(nimble)
library(ggmcmc)
library(parallel)
library(MCMCvis)

### data: weekly_data_Cuba_COVID_2021_population_lag

ID_municipalities <- weekly_data_Cuba_COVID_2021_population_lag$ID
Population_size_municipalities <- weekly_data_Cuba_COVID_2021_population_lag$Total
Timepoints_weeks <- weekly_data_Cuba_COVID_2021_population_lag$week
y_deaths <- weekly_data_Cuba_COVID_2021_population_lag$deaths_lag2
y_imported_cases <- weekly_data_Cuba_COVID_2021_population_lag$Imported

e = structure(.Data = Population_size_municipalities, .Dim = c(44, 168))
e = t(e)
y_d = structure(.Data = y_deaths, .Dim = c(44, 168))
y_d = t(y_d)
y_ic = structure(.Data = y_imported_cases, .Dim = c(44, 168))
y_ic = t(y_ic)

T <- 44
m <- 168

nbWB_A <- nb2WB(nb = nb_q)
adj <- nbWB_A$adj
adj
num <- nbWB_A$num
sumNumNeigh <-  sum(num)
sumNumNeigh
L <- sum(num)

R = structure(.Data = c(1, 0, 0, 1), .Dim = c(2,2))
R

Q = structure(.Data = c(1, 0, 0, 1), .Dim = c(2,2))
Q

mu_g = rep(0, 2)
mu_g

#constants
LSTConsts_3 <-list(T, m, adj, sumNumNeigh, num, L, e, R, Q, mu_g)
names(LSTConsts_3) <-c("T","m","adj","sumNumNeigh","num", "L", "e", "R", "Q", "mu_g")


#data
LSTData <- list(y_d, y_ic)
names(LSTData) <- c("y_d", "y_ic")


############################################
############## NIMBLE MODEL  ###############
############################################
detectCores()
nbcores <- 4
cluster1 <- makeCluster(nbcores)
workflow_model <- function(seed, data, cst, inits, prms) {
  
  library(nimble)
  
  STCode_model <-nimbleCode(
    {
      for (i in 1:m)
      {
        for (k in 1:T)
        { y_d[i,k]~dpois(mu_d[i,k])
          pred_y_d[i,k]~dpois(mu_d[i,k])
          log(mu_d[i,k])<-log(e[i,k])+log(theta_d[i,k])
          log(theta_d[i,k])<-a0_d+g[k, 1]+phi[i,1]+v_d[i]+psi_d[i,k]
          psi_d[i,k]~dnorm(0,taupsi_d)
          Dp_d[i,k]<-(y_d[i,k]- pred_y_d[i,k])**2
          
          y_ic[i,k]~dpois(mu_ic[i,k])
          pred_y_ic[i,k]~dpois(mu_ic[i,k])
          log(mu_ic[i,k])<-log(e[i,k])+log(theta_ic[i,k])
          log(theta_ic[i,k])<-a0_ic+g[k, 2]+phi[i,2]+v_ic[i]+psi_ic[i,k]
          psi_ic[i,k]~dnorm(0,taupsi_ic)
          Dp_ic[i,k]<-(y_ic[i,k]- pred_y_ic[i,k])**2
          
          exp_g[k, 1] <- exp(g[k, 1])
          exp_g[k, 2] <- exp(g[k, 2])
          
          exp_phi[i,1] <- exp(phi[i,1]) 
          exp_phi[i,2] <- exp(phi[i,2]) 
        }	
        v_d[i]~dnorm(0,tauv_d)
        v_ic[i]~dnorm(0,tauv_ic)
      }
      
      g[1, 1:2]~dmnorm(mu_g[1:2],inv_D[1:2, 1:2])
      
      for (k in 2:T){
        for (j in 1:2){
          mu_vector_g[k,j] <- g[k-1,j]
        }
        g[k, 1:2]~dmnorm(mu_vector_g[k,1:2],inv_D[1:2, 1:2])
      }
      
      inv_D[1:2, 1:2] ~ dwish(Q[1:2,1:2], 2) 
      sigma2.g[1:2, 1:2] <- inverse(inv_D[1:2,1:2]) 
      
      sigma.g[1] <- sqrt(sigma2.g[1, 1])
      sigma.g[2] <- sqrt(sigma2.g[2, 2])
      corr.g <- sigma2.g[1, 2] / (sigma.g[1] * sigma.g[2]) 
      
      
      for (k in 1:2)
      {
        u[1:m,k]~dcar_normal(adj[1:L],weights[1:L],num[1:m],tauu[k],zero_mean=1)	
      }
      
      for(k in 1:sumNumNeigh) {
        weights[k] <- 1
      }
      
      Prec[1:2,1:2]~dwish(R[1:2,1:2],2)
      Cov[1:2,1:2]<-inverse(Prec[1:2,1:2])
      Achol[1:2,1:2]<-chol(Cov[1:2,1:2])
      
      for(i in 1:m) {
        phi[i,1:2]<-Achol[1:2,1:2]%*%u[i,1:2]
      }
      
      tauu[1]<-1/pow(sdu_d,2)
      sdu_d~dunif(0,10)
      tauu[2]<-1/pow(sdu_ic,2)
      sdu_ic~dunif(0,10)
      
      tauv_d<-1/pow(sdv_d,2)
      sdv_d~dunif(0,10)
      tauv_ic<-1/pow(sdv_ic,2)
      sdv_ic~dunif(0,10)
      
      taupsi_d<-1/pow(sdpsi_d,2)
      sdpsi_d~dunif(0,10)
      taupsi_ic<-1/pow(sdpsi_ic,2)
      sdpsi_ic~dunif(0,10)
      
      a0_d~dflat()
      a0_ic~dflat()
      
      sig[1]<-sqrt(Cov[1,1])
      sig[2]<-sqrt(Cov[2,2])
      
      cor12<-Cov[1,2]/(sig[1]*sig[2])    
      
      mspe_d <-mean(Dp_d[1:m,1:T])
      mspe_ic <-mean(Dp_ic[1:m,1:T])
      
    })
  
  RW1_int <- nimbleModel(code = STCode_model, 
                             data = data, 
                             inits = inits,
                             constants = cst)
  C_RW1_int <- compileNimble(RW1_int)
  RW1_intConf <- configureMCMC(RW1_int)
  RW1_intConf$addMonitors(c("exp_phi", "exp_g", "theta_d", "theta_ic",
                                "corr.g", "cor12", "sigma2.g", "sigma.g",
                                "Cov", "sig", "mspe_d", "mspe_ic"))
  RW1_intMCMC <- buildMCMC(RW1_intConf, project=RW1_int)
  
  C_RW1_intMCMC <- compileNimble(RW1_intMCMC)
  
  RW1_int_samples <- runMCMC(mcmc = C_RW1_intMCMC, 
                                 niter = 1100000, 
                                 nburnin = 100000,
                                 thin = 1000,
                                 setSeed = seed)
  
  
  return(RW1_int_samples)
}

LSTinits_model <-list(a0_d=-12,
                           sdv_d=0.1,
                           sdu_d=0.1, sdpsi_d=0.1,
                           v_d=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                           psi_d=structure(.Data=rep(0.1,7392),.Dim=c(168,44)),
                           pred_y_d=structure(.Data=rep(0,7392),.Dim=c(168,44)),
                           
                           a0_ic=-12,
                           sdv_ic=0.1,
                           sdpsi_ic=0.1, sdu_ic = 0.1,
                           v_ic=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                           u = structure(.Data=rep(0,336),.Dim=c(168,2)),
                           g = structure(.Data=rep(0,88),.Dim=c(44,2)),
                           pred_y_ic=structure(.Data=rep(0,7392),.Dim=c(168,44)),
                           psi_ic=structure(.Data=rep(0.1,7392),.Dim=c(168,44)),
                           Prec=structure(.Data=c(1, 0,0,1),.Dim=c(2,2)),
                           Inv_D=structure(.Data=c(1, 0,0,1),.Dim=c(2,2)))

params_model <- c("a0_d", "sdu_d", "sdv_d", "sdpsi_d",
                       "a0_ic", "sdv_ic", "sdpsi_ic", "sdu_ic",
                       "cor12",
                       "corr.g", "sigma2.g")

output_model <- parLapply(cl = cluster1, 
                               X = c(222, 666, 444, 555),
                               fun = workflow_model,
                               data = LSTData,
                               cst = LSTConsts,
                               inits = LSTinits_model,
                               prms = params_model)
output_model
MCMCsummary(output_model, params = c("corr.g", "cor12", "a0_d", "a0_ic", 
                                          "sdu_ic", "sdu_d",
                                          "sdpsi_ic", "sdpsi_d",
                                          "sdv_ic", "sdv_d", "corr.g", "cor12", 
                                          "sigma2.g", "sigma.g",
                                          "Cov", "sig", "mspe_d", "mspe_ic"))

MCMCtrace(output_model, params = c("corr.g", "cor12", "a0_d", "a0_ic", 
                                          "sdu_ic", "sdu_d",
                                          "sdpsi_ic", "sdpsi_d",
                                          "sdv_ic", "sdv_d"))
stopCluster(cluster1)

#PLOT TEMPORAL TREND G
#DEATHS
MCMCsummary(output_model, params = "exp_g")
matrix_summary <- as.matrix(MCMCsummary(output_model, params = "exp_g"))
temporal_means <- matrix_summary[1:44,c(1,3,5)]
temporal_means
weeks_g <- 1:44
weeks_g
temporal_means_df <- as.data.frame(temporal_means)
temporal_means_df$weeks <- weeks_g
temporal_means_df
colnames(temporal_means_df)[2] <- "Lower"
colnames(temporal_means_df)[3] <- "Upper"

ggplot(data=temporal_means_df, aes(x=weeks, y=mean, group=1)) +
  geom_line()

ggplot(temporal_means_df, aes(x=weeks, group=1)) + 
  geom_line(aes(y = mean), color = "black") + 
  geom_line(aes(y = Lower), color="black", linetype="twodash") + 
  geom_line(aes(y = Upper), color="black", linetype="twodash") + 
  labs(x = "Week", y = "exp(g)", title = "Overall temporal trend of deaths (lag 6)")

#IMPORTED CASES
MCMCsummary(output_model, params = "exp_g")
matrix_summary <- as.matrix(MCMCsummary(output_model, params = "exp_g"))
temporal_means <- matrix_summary[45:88,c(1,3,5)]
temporal_means
weeks_g <- 1:44
weeks_g
temporal_means_df <- as.data.frame(temporal_means)
temporal_means_df$weeks <- weeks_g
temporal_means_df
colnames(temporal_means_df)[2] <- "Lower"
colnames(temporal_means_df)[3] <- "Upper"

ggplot(data=temporal_means_df, aes(x=weeks, y=mean, group=1)) +
  geom_line()

ggplot(temporal_means_df, aes(x=weeks, group=1)) + 
  geom_line(aes(y = mean), color = "black") + 
  geom_line(aes(y = Lower), color="black", linetype="twodash") + 
  geom_line(aes(y = Upper), color="black", linetype="twodash") + 
  labs(x = "Week", y = "exp(g)", title = "Overall temporal trend of imported cases")

#PLOT SPATIAL TREND U
#DEATHS
MCMCsummary(output_model, params = "exp_phi")
matrix_summary <- as.matrix(MCMCsummary(output_model, params = "exp_phi"))
matrix_summary
spatial_re <- matrix_summary[1:168,1]
spatial_re
municipalities_id <- 1:168
municipalities_id
spatial_re_df <- as.data.frame(spatial_re)
spatial_re_df$ID <- municipalities_id
spatial_re_df

cuba_df <- cuba_shape %>%
  left_join(spatial_re_df, by = c("ID" = "ID"))

ggplot(data = cuba_df) +
  geom_sf(aes(fill = spatial_re), show.legend = NA)+
  scale_fill_gradient(name="exp(v)")+ 
  labs(title = "Spatial trend of deaths")

#IMPORTED CASES
MCMCsummary(output_model, params = "exp_phi")
matrix_summary <- as.matrix(MCMCsummary(output_model, params = "exp_phi"))
matrix_summary
spatial_re <- matrix_summary[169:336,1]
spatial_re
municipalities_id <- 1:168
municipalities_id
spatial_re_df <- as.data.frame(spatial_re)
spatial_re_df$ID <- municipalities_id
spatial_re_df

cuba_df <- cuba_shape %>%
  left_join(spatial_re_df, by = c("ID" = "ID"))

ggplot(data = cuba_df) +
  geom_sf(aes(fill = spatial_re), show.legend = NA)+
  scale_fill_gradient(name="exp(v)")+ 
  labs(title = "Spatial trend of imported cases")

### CALCULATE WAIC ###
#global
STCode_model <-nimbleCode(
  {
    for (i in 1:m)
    {
      for (k in 1:T)
      { y_d[i,k]~dpois(mu_d[i,k])
        pred_y_d[i,k]~dpois(mu_d[i,k])
        log(mu_d[i,k])<-log(e[i,k])+log(theta_d[i,k])
        log(theta_d[i,k])<-a0_d+g[k, 1]+phi[i,1]+v_d[i]+psi_d[i,k]
        psi_d[i,k]~dnorm(0,taupsi_d)
        Dp_d[i,k]<-(y_d[i,k]- pred_y_d[i,k])**2
        
        y_ic[i,k]~dpois(mu_ic[i,k])
        pred_y_ic[i,k]~dpois(mu_ic[i,k])
        log(mu_ic[i,k])<-log(e[i,k])+log(theta_ic[i,k])
        log(theta_ic[i,k])<-a0_ic+g[k, 2]+phi[i,2]+v_ic[i]+psi_ic[i,k]
        psi_ic[i,k]~dnorm(0,taupsi_ic)
        Dp_ic[i,k]<-(y_ic[i,k]- pred_y_ic[i,k])**2
        
        exp_g[k, 1] <- exp(g[k, 1])
        exp_g[k, 2] <- exp(g[k, 2])
        
        exp_phi[i,1] <- exp(phi[i,1]) 
        exp_phi[i,2] <- exp(phi[i,2]) 
      }	
      v_d[i]~dnorm(0,tauv_d)
      v_ic[i]~dnorm(0,tauv_ic)
    }
    
    g[1, 1:2]~dmnorm(mu_g[1:2],inv_D[1:2, 1:2])
    
    for (k in 2:T){
      for (j in 1:2){
        mu_vector_g[k,j] <- g[k-1,j]
      }
      g[k, 1:2]~dmnorm(mu_vector_g[k,1:2],inv_D[1:2, 1:2])
    }
    
    inv_D[1:2, 1:2] ~ dwish(Q[1:2,1:2], 2) 
    sigma2.g[1:2, 1:2] <- inverse(inv_D[1:2,1:2]) 
    
    sigma.g[1] <- sqrt(sigma2.g[1, 1])
    sigma.g[2] <- sqrt(sigma2.g[2, 2])
    corr.g <- sigma2.g[1, 2] / (sigma.g[1] * sigma.g[2]) 
    
    
    for (k in 1:2)
    {
      u[1:m,k]~dcar_normal(adj[1:L],weights[1:L],num[1:m],tauu[k],zero_mean=1)	
    }
    
    for(k in 1:sumNumNeigh) {
      weights[k] <- 1
    }
    
    Prec[1:2,1:2]~dwish(R[1:2,1:2],2)
    Cov[1:2,1:2]<-inverse(Prec[1:2,1:2])
    Achol[1:2,1:2]<-chol(Cov[1:2,1:2])
    
    for(i in 1:m) {
      phi[i,1:2]<-Achol[1:2,1:2]%*%u[i,1:2]
    }
    
    tauu[1]<-1/pow(sdu_d,2)
    sdu_d~dunif(0,10)
    tauu[2]<-1/pow(sdu_ic,2)
    sdu_ic~dunif(0,10)
    
    tauv_d<-1/pow(sdv_d,2)
    sdv_d~dunif(0,10)
    tauv_ic<-1/pow(sdv_ic,2)
    sdv_ic~dunif(0,10)
    
    taupsi_d<-1/pow(sdpsi_d,2)
    sdpsi_d~dunif(0,10)
    taupsi_ic<-1/pow(sdpsi_ic,2)
    sdpsi_ic~dunif(0,10)
    
    a0_d~dflat()
    a0_ic~dflat()
    
    sig[1]<-sqrt(Cov[1,1])
    sig[2]<-sqrt(Cov[2,2])
    
    cor12<-Cov[1,2]/(sig[1]*sig[2])    # between deaths and imported cases
    
    mspe_d <-mean(Dp_d[1:m,1:T])
    mspe_ic <-mean(Dp_ic[1:m,1:T])
    
  })

RW1_int <- nimbleModel(code = STCode_model, 
                           data = LSTData, 
                           inits = LSTinits_model,
                           constants = LSTConsts)

C_RW1_int <- compileNimble(RW1_int)
RW1_intConf <- configureMCMC(RW1_int)
RW1_intConf$addMonitors(c("exp_phi", "exp_g", "theta_d", "theta_ic",
                              "corr.g", "cor12", "sigma2.g", "sigma.g",
                              "Cov", "sig", "mspe_d", "mspe_ic",
                              "u", "v_d", "v_ic", "psi_d", "psi_ic", "g"))
RW1_intMCMC <- buildMCMC(RW1_intConf, project=RW1_int)

C_RW1_intMCMC <- compileNimble(RW1_intMCMC)

posteriorSamplesMatrix <- rbind(output_model[[1]], output_model[[2]], 
                                output_model[[3]], output_model[[4]])
C_RW1_intMCMC$run(5)   
nimble:::matrix2mv(posteriorSamplesMatrix, C_RW1_intMCMC$mvSamples)
C_RW1_intMCMC$enableWAIC <- TRUE
calculateWAIC(C_RW1_intMCMC)