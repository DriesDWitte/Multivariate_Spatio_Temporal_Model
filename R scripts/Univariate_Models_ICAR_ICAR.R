##############################################################################
# This script contains the code to fit the univariate spatio-temporal model  #
# using the ICAR prior for the spatial random effects and the ICAR prior     #
# for temporal random effects.                                               #
##############################################################################

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
y_deaths <- weekly_data_Cuba_COVID_2021_population_lag$Deaths
y_imported_cases <- weekly_data_Cuba_COVID_2021_population_lag$Imported

e = structure(.Data = Population_size_municipalities, .Dim = c(44, 168))
e = t(e)
y = structure(.Data = y_deaths, .Dim = c(44, 168))
y = t(y)

T <- 44
m <- 168

nbWB_A <- nb2WB(nb = nb_q)
adj <- nbWB_A$adj
adj
num <- nbWB_A$num
sumNumNeigh <-  sum(num)
sumNumNeigh
L <- sum(num)

Wg <- matrix(0, nrow = 44, ncol = 44)

for(i in 1:(44-1)) Wg[i,i+1] <- 1
for(i in 1:(44-1)) Wg[i+1,i] <- 1

Wgnb <- mat2listw(Wg)
Wgnb <- nb2WB(nb = Wgnb$neighbours)

L_g = length(Wgnb$adj)
adj_g = Wgnb$adj                        
num_g = Wgnb$num
weights_g = Wgnb$weights

pop_density <- weekly_data_Cuba_COVID_2021_population_lag$pop_density_devided

x = structure(.Data = pop_density, .Dim = c(44, 168))
x = t(x)

#constants
LSTConsts <-list(T, m, adj, sumNumNeigh, num, L, e, L_g, adj_g, num_g, weights_g, x)
names(LSTConsts) <-c("T","m","adj","sumNumNeigh","num", "L", "e", "L_g",
                       "adj_g", "num_g", "weights_g", "x")
LSTConsts


#data
y = structure(.Data = y_deaths, .Dim = c(44, 168))
y = t(y)
LSTData <- list(y)
names(LSTData) <- c("y")


############################################
############## NIMBLE MODEL  ###############
############################################

#make cluster
detectCores()
nbcores <- 4
my_cluster <- makeCluster(nbcores)

LSTinits_RW1_int <-list(a0=-12, b1=0,
                            sdv=0.1,
                            sdg=0.1,sdu=0.1, sdpsi=0.1,
                            g=c(0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0),
                            v=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                            u=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                            pred_y=structure(.Data=rep(0,7392),.Dim=c(168,44)),
                            psi=structure(.Data=rep(0.1,7392),.Dim=c(168,44)))

params_RW1_int <- c("a0", "sdu", "sdv", "sdg", "sdpsi",
                        "exp_u", "exp_g", "u_exced", "b1", "mspe")

workflow_RW1_int <- function(seed, data, cst, inits, prms) {
  
  library(nimble)
  
  RW1_int_model <-nimbleCode(
    {
      for (i in 1:m)
      {
        for (k in 1:T)
        {y[i,k]~dpois(mu[i,k])
          pred_y[i,k]~dpois(mu[i,k])
          log(mu[i,k])<-log(e[i,k])+log(theta[i,k])
          log(theta[i,k])<-a0+g[k]+u[i]+v[i]+psi[i,k]+b1*x[i,k]
          psi[i,k]~dnorm(0,taupsi)
          Dp[i,k]<-(y[i,k]- pred_y[i,k])**2
          
          exp_g[k] <- exp(g[k]) 
         
        }	
        v[i]~dnorm(0,tauv)
        
        exp_u[i] <- exp(u[i]) 
        exp_v[i] <- exp(v[i]) 
        
        #exceedance probabilities
        u_exced[i]<-step(exp_u[i]-1)
        v_exced[i]<-step(exp_v[i]-1)
      }
      
      g[1:T]~dcar_normal(adj_g[1:L_g],weights_g[1:L_g],num_g[1:T],taug, zero_mean=1)
      
      u[1:m]~dcar_normal(adj[1:L],weights[1:L],num[1:m],tauu, zero_mean=1)
      
      for(k in 1:sumNumNeigh) {
        weights[k] <- 1
      }
      
      tauu<-1/pow(sdu,2)
      sdu~dunif(0,10)
      taug<-1/pow(sdg,2)
      sdg~dunif(0,10)
      tauv<-1/pow(sdv,2)
      sdv~dunif(0,10)
      taupsi<-1/pow(sdpsi,2)
      sdpsi~dunif(0,10)
      
      a0~dflat()
      b1~dnorm(0, sd = 100)
      
      mspe <-mean(Dp[1:m,1:T])
      
    })
  

  
  RW1_int <- nimbleModel(code = RW1_int_model, 
                    data = data, 
                    inits = inits,
                    constants = cst)
  C_RW1_int <- compileNimble(RW1_int)
  RW1_intConf <- configureMCMC(RW1_int)
  RW1_intConf$addMonitors(c("exp_u", "exp_g", "exp_v", "theta", "mspe", 
                                "u_exced", "v_exced"))
  RW1_intMCMC <- buildMCMC(RW1_intConf, project=RW1_int)
  
  C_RW1_intMCMC <- compileNimble(RW1_intMCMC)
  
  RW1_int_samples <- runMCMC(mcmc = C_RW1_intMCMC, 
                            niter = 110000, 
                            nburnin = 10000,
                            thin = 100,
                            setSeed = seed)
  
  return(RW1_int_samples)
}

output <- parLapply(cl = my_cluster, 
                    X = c(222, 666, 444, 555),
                    fun = workflow_RW1_int,
                    data = LSTData,
                    cst = LSTConsts,
                    inits = LSTinits_RW1_int,
                    prms = params_RW1_int)
output
options(max.print=10000)
options(scipen=999)
MCMCsummary(output, params = c("a0", "sdg", "sdv", "sdu", "sdpsi", "mspe", "b1"))
MCMCsummary(output, params = "b1")
MCMCtrace(output, params = "b1")

stopCluster(my_cluster)

#PLOT TEMPORAL TREND G
#obtain mean of exp_g
matrix_summary <- as.matrix(MCMCsummary(output, params = "exp_g"))
temporal_means <- matrix_summary[,c(1,3,5)]
temporal_means
weeks_g <- 1:44
weeks_g
temporal_means_df <- as.data.frame(temporal_means)
temporal_means_df$weeks <- weeks_g
temporal_means_df
colnames(temporal_means_df)[2] <- "Lower"
colnames(temporal_means_df)[3] <- "Upper"

ggplot(temporal_means_df, aes(x=weeks, group=1)) + 
  geom_line(aes(y = mean), color = "black") + 
  geom_line(aes(y = Lower), color="black", linetype="twodash") + 
  geom_line(aes(y = Upper), color="black", linetype="twodash") + 
  labs(x = "Week", y = "exp(g)", title = "Temporal trend of the number of deaths")

#PLOT SPATIAL TREND U
#obtain mean of exp_u
matrix_summary <- as.matrix(MCMCsummary(output, params = "exp_u"))
matrix_summary
spatial_re <- matrix_summary[,1]
spatial_re
municipalities_id <- 1:168
municipalities_id
spatial_re_df <- as.data.frame(spatial_re)
spatial_re_df$ID <- municipalities_id
spatial_re_df

cuba_df <- cuba_shape %>%
  left_join(spatial_re_df, by = c("ID" = "ID"))

data <- cuba_df[with(cuba_df,order(-spatial_re)),]

data <- data[1:10,]
data

ggplot(data = cuba_df) +
  geom_sf(aes(fill = spatial_re), show.legend = NA)+
  scale_fill_gradient(name="exp(u)", low = "white", high = "gray0")+ 
  labs(title = "Spatial trend of the number deaths")


#obtain mean of exced_u
matrix_summary <- as.matrix(MCMCsummary(output, params = "u_exced"))
matrix_summary
spatial_re <- matrix_summary[,1]
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
  scale_fill_gradient(name="P(exp(u) > 1)", low = "white", high = "gray0")+ 
  labs(title = "Exceedance probabilities for the number of deaths")

test <- cuba_df %>% mutate(new_bin = cut(spatial_re, breaks=c(0, 0.05, 0.95, 1)))

ggplot(data = test) +
  geom_sf(aes(fill = factor(new_bin)), show.legend = NA)+
  #scale_fill_gradient(name="exp(v)")+ 
  theme_bw()+
  scale_fill_manual(values = c("grey100", "grey50", "grey0"), name= "exp(v)")+ 
  labs(title = "Exceedence probabilities of Deaths (structured)")

#PLOT SPATIAL TREND V
#obtain mean of exp_v
matrix_summary <- as.matrix(MCMCsummary(output, params = "exp_v"))
matrix_summary
spatial_re <- matrix_summary[,1]
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
  scale_fill_gradient(name="exp(u)", low = "white", high = "gray0")+ 
  labs(title = "Overall spatial trend of deaths (unstructured)")

#obtain mean of exced_v
matrix_summary <- as.matrix(MCMCsummary(output, params = "v_exced"))
matrix_summary
spatial_re <- matrix_summary[,1]
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
  scale_fill_gradient(name="exp(v)", low = "white", high = "gray0")+ 
  labs(title = "Exceedence probabilities of Deaths (unstructured)")

test <- cuba_df %>% mutate(new_bin = cut(spatial_re, breaks=c(0, 0.05, 0.95, 1)))

ggplot(data = test) +
  geom_sf(aes(fill = factor(new_bin)), show.legend = NA)+
  #scale_fill_gradient(name="exp(v)")+ 
  theme_bw()+
  scale_fill_manual(values = c("grey100", "grey50", "grey0"), name= "exp(v)")+ 
  labs(title = "Exceedence probabilities of Deaths (unstructured)")

### CALCULATE WAIC ###
#global
RW1_int_model <-nimbleCode(
    {
      for (i in 1:m)
      {
        for (k in 1:T)
        {y[i,k]~dpois(mu[i,k])
          pred_y[i,k]~dpois(mu[i,k])
          log(mu[i,k])<-log(e[i,k])+log(theta[i,k])
          log(theta[i,k])<-a0+g[k]+u[i]+v[i]+psi[i,k]
          psi[i,k]~dnorm(0,taupsi)
          Dp[i,k]<-(y[i,k]- pred_y[i,k])**2
          
          exp_g[k] <- exp(g[k]) 
          
        }	
        v[i]~dnorm(0,tauv)
        
        exp_u[i] <- exp(u[i]) 
        exp_v[i] <- exp(v[i]) 
        
        #exceedance probabilities
        u_exced[i]<-step(exp_u[i]-1)
        v_exced[i]<-step(exp_v[i]-1)
      }
      
      g[1:T]~dcar_normal(adj_g[1:L_g],weights_g[1:L_g],num_g[1:T],taug, zero_mean=1)
      
      u[1:m]~dcar_normal(adj[1:L],weights[1:L],num[1:m],tauu, zero_mean=1)
      
      for(k in 1:sumNumNeigh) {
        weights[k] <- 1
      }
      
      tauu<-1/pow(sdu,2)
      sdu~dunif(0,10)
      taug<-1/pow(sdg,2)
      sdg~dunif(0,10)
      tauv<-1/pow(sdv,2)
      sdv~dunif(0,10)
      taupsi<-1/pow(sdpsi,2)
      sdpsi~dunif(0,10)
      
      a0~dflat()
      mspe <-mean(Dp[1:m,1:T])
      
    })
  
  RW1_int <- nimbleModel(code = RW1_int_model, 
                             data = LSTData, 
                             inits = LSTinits_RW1_int,
                             constants = LSTConsts)
  C_RW1_int <- compileNimble(RW1_int)
  RW1_intConf <- configureMCMC(RW1_int)
  RW1_intConf$addMonitors(c("exp_u", "exp_g", "exp_v", "theta", "mspe", 
                                "u_exced", "v_exced", "u", "g", "v", "psi"))
  RW1_intMCMC <- buildMCMC(RW1_intConf, project=RW1_int)
  
  C_RW1_intMCMC <- compileNimble(RW1_intMCMC)
  

posteriorSamplesMatrix <- rbind(output[[1]], output[[2]], output[[3]], output[[4]])
C_RW1_intMCMC$run(5)   ## non-zero number of iterations
nimble:::matrix2mv(posteriorSamplesMatrix, C_RW1_intMCMC$mvSamples)
C_RW1_intMCMC$enableWAIC <- TRUE
calculateWAIC(C_RW1_intMCMC)
