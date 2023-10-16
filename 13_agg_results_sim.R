###
### Aggregating simulation results
###

# rm(list=ls())

setwd("~/Documents/ipm/cwd_ipm_sim4_chtc")

library(tidyverse)
library(nimble)
library(Matrix)
library(coda)
library(lattice)


####################################################
###
### Loading True Values
###
####################################################

# source("02_load_all_data_to_run.R")

###############################################################
# Generate simulated data
###############################################################

#setting seed: first processnum == 0
# processnum <- 100
# seedy <- processnum + 1000
# set.seed(seedy)

# source("04_generate_data.R")

###
### read in results from CHTC runs
###

### number of simulations run on CHTC
n_sims <- 100

### initiating vectors to save output 
### from each simulation process

fit_sum_ls <- vector(mode = "list", length = n_sims)
gd_ls <- vector(mode = "list", length = n_sims)
runtime_ls <- vector(mode = "list", length = n_sims)
seed_ls <- vector(mode = "list", length = n_sims)
ess_ls <- vector(mode = "list", length = n_sims)
# true.values.ls=vector(mode = "list", length = n_sims)


modelid <- "A"
sims <- 1:n_sims
#if not done running all
# sims <- (1:n_sims)[-c(13)]

for(i in sims){
  load(paste0('results/',modelid,'/fit_sum_',modelid,'_',i,'.Rdata'))
  fit_sum_ls[[i]] <- fit_sum

  load(paste0('results/',modelid,'/gd_',modelid,'_',i,'.Rdata'))
  gd_ls[[i]] <- gd

  load(paste0('results/',modelid,'/runtime_',modelid,'_',i,'.Rdata'))
  runtime_ls[[i]] <- runtime

  load(paste0('results/',modelid,'/seed_',modelid,'_',i,'.Rdata'))
  seed_ls[[i]] <- seedy

  load(paste0('results/',modelid,'/ess_',modelid,'_',i,'.Rdata'))
  ess_ls[[i]] <- ess
}

###
### combine values across lists
###

runtime <- unlist(runtime_ls)

seeds <- unlist(seed_ls)

runtime_mn <- mean(runtime,na.rm=TRUE)

# ess.df  <- data.frame(do.call(rbind,ess_ls))
# names(ess.df)
# names(ess.df)="beta0"
# ess.df$simulation=1:n_sims

################################################
# true values to calculate MSE and Coverage,
# these are the same across all simulations
################################################

beta0_survival_inf_true
beta0_survival_sus_true
nT_age
nT_period
age_effect_true
period_effect_true
### this will change in the next run
foi_age_true <- m_age_foi
n_ageclass <- n_ageclassm
sn_sus_true <- c(sn_sus_true)
sn_inf_true <- c(sn_inf_true)

######################################################
###
### Susceptile intercept Beta0
### Calculating Bias, Mean Sq Error, and Coverage
### 
######################################################

beta0_survival_sus_bias<-rep(NA,n_sims)
beta0_survival_sus_coverage_ind<-rep(NA,n_sims)

for(i in sims){
  temp <- fit_sum_ls[[i]][grep("beta0_survival_sus",rownames(fit_sum_ls[[i]])),]
  beta0_survival_sus_mn <- temp[1]
  beta0_survival_sus_bias[i] <- beta0_survival_sus_mn - beta0_survival_sus_true
  beta0_survival_sus_coverage_ind[i] <- ifelse(beta0_survival_sus_true >temp[4] & beta0_survival_sus_true <temp[5],1,0)
}

beta0_survival_sus_mse <- (1/length(sims)) * sum(beta0_survival_sus_bias^2,na.rm = TRUE)

beta0_survival_sus_mse

png(paste0("figures/",modelid,"/beta0_survival_sus_bias_hist_",modelid,".png"))
  hist(beta0_survival_sus_bias, breaks = 30)
  abline(v = 0,lty = 2,lwd=3)
dev.off()

beta0_survival_sus_coverage <- sum(beta0_survival_sus_coverage_ind,na.rm = TRUE)/length(sims)



######################################################
###
### Infected intercept Beta0
### Calculating Bias, Mean Sq Error, and Coverage
### 
######################################################

beta0_survival_inf_bias <- rep(NA, n_sims)
beta0_survival_inf_coverage_ind <- rep(NA, n_sims)

for(i in sims){
  temp <- fit_sum_ls[[i]][grep("beta0_survival_inf",rownames(fit_sum_ls[[i]])),]
  beta0_survival_inf_mn <- temp[1]
  beta0_survival_inf_bias[i] <- beta0_survival_inf_mn - beta0_survival_inf_true
  beta0_survival_inf_coverage_ind[i] <- ifelse(beta0_survival_inf_true >temp[4] & beta0_survival_inf_true <temp[5],1,0)
}

beta0_survival_inf_mse <- (1/length(sims)) * sum(beta0_survival_inf_bias^2, na.rm = TRUE)


png(paste0("figures/",modelid,"/beta0_survival_inf_bias_hist_",modelid,".png"))
  hist(beta0_survival_inf_bias, breaks = 30)
  abline(v = 0,lty = 2,lwd=3)
dev.off()
beta0_survival_inf_coverage <- sum(beta0_survival_inf_coverage_ind,na.rm = TRUE)/n_sims


beta0_tab <- round(data.frame(beta0_survival_sus = c(beta0_survival_sus_mse,beta0_survival_sus_coverage),
  beta0_survival_inf = c(beta0_survival_inf_mse,beta0_survival_inf_coverage)),3)

rownames(beta0_tab) <- c("MSE","Coverage")

write.csv(beta0_tab, file =paste0("results/",modelid,"/beta0_survival_tab_",modelid,".csv"))

#########################################################
### Color pallets
#########################################################
library(MetBrewer)
renoir_pal <- met.brewer(name="Renoir", n=10, type="discrete")
pillement_pal <- met.brewer(name="Pillement", n=6, type="discrete")
troy_pal <- met.brewer(name="Troy", n=8, type="discrete")


###
### plot of bias density
###

df_beta0_bias  <- data.frame(Susceptible = beta0_survival_sus_bias,Infected = beta0_survival_inf_bias) %>% pivot_longer(cols = everything())
df_beta0_bias$name <- as.factor(df_beta0_bias$name)
df_beta0_bias$name <- factor(df_beta0_bias$name,levels = c('Susceptible','Infected'))




beta0_bias_density_plot <- ggplot(data = df_beta0_bias) +
  geom_density(aes(x = value,color = name,fill = name),alpha = .9) +
  geom_vline(xintercept = 0, linetype = "dashed",size = 1.5) +
  theme_bw()+
  ggtitle("Bias of Mortality Hazard Intercepts")+
  xlab("Bias")+
  ylab("Density")+
  xlim(-1.2,0.6)+
  scale_fill_manual("CWD Status",values = troy_pal[c(2, 7)])+
  scale_color_manual("CWD Status",values = troy_pal[c(2, 7)])+
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )

beta0_bias_density_plot

ggsave(paste0("figures/",modelid,"/beta0_bias_density_",modelid,".png"),
      beta0_bias_density_plot,
      height = 6,
      width = 10)



#####################################################
###
### Age effects bias and coverage
###
#####################################################

age_effect_bias <- vector(mode = "list", length = n_sims)
age_effect_coverage_ind <- vector(mode = "list", length = n_sims)

#removing the couple of seeds that failed for now
#fixing data generating should fix this too
# sims <- 1:n_sims
# sims <- (1:n_sims)[-3]
# sims <- (sims)[-29]


for(i in sims) {
  temp <- fit_sum_ls[[i]][grep("age_effect_survival",rownames(fit_sum_ls[[i]])),]
  age_effect_mn <- temp[, 1]
  age_effect_bias[[i]] <- age_effect_mn - age_effect_true
  coverage_ind <- rep(NA, nT_age)
  for(j in 1:nT_age){
    coverage_ind[j] <- ifelse(age_effect_true[j]>temp[j,4] & age_effect_true[j] < temp[j,5],1,0)
  }
  age_effect_coverage_ind[[i]] <- coverage_ind
}

age_effect_coverage_mat <- do.call(rbind, age_effect_coverage_ind)
age_indicators <- apply(age_effect_coverage_mat, 2, sum)
age_effect_coverage <- age_indicators / length(sims)
age_effect_bias_df <- do.call(rbind, age_effect_bias)
age_effect_mse <- apply(age_effect_bias_df^2, 2, mean)


plot(age_effect_mse)
plot(age_effect_coverage)

##################################
### Proper plots of age effects
##################################

df_age_effect_out <- data.frame(mse = age_effect_mse,
                                  coverage = age_effect_coverage)

df_age_effect_out$age <- 1:nT_age
df_age_effect_out <- df_age_effect_out %>% pivot_longer(cols = c('mse','coverage')) 
df_age_effect_out$name <- as.factor(df_age_effect_out$name)
levels(df_age_effect_out$name) <- c("Coverage","MSE")

#########################
### age effect mse
#########################
age_mse_coverage_plot <- ggplot(data = df_age_effect_out)+
  geom_line(aes(x = age,y = value,color = name),alpha = .9, size = 1.5)+
  facet_wrap(.~name,scales = "free") +
  theme_bw()+
  ggtitle("Age Effects Results")+
  ylab("")+
  xlab("Age(Years)")+
  scale_color_manual(values = renoir_pal[c(1, 6)])+
  scale_x_continuous(breaks = seq(104,nT_age,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
age_mse_coverage_plot

ggsave(paste0("figures/",modelid,"/age_mse_coverage_plot_",modelid,".png"),
      age_mse_coverage_plot,
      height = 6,
      width = 8)

df_age_effect_long <-  data.frame(mse = age_effect_mse,
                                  coverage = age_effect_coverage) %>% pivot_longer(cols = everything())
df_age_effect_long$name <- as.factor(df_age_effect_long$name)
levels(df_age_effect_long$name) <- c("Coverage","MSE")

age_effect_mse_plot <- ggplot(data = df_age_effect_long)+
  geom_density(aes(x = value,color = name,fill = name),alpha = .9)+
  facet_wrap(.~name,scales = "free") +
  theme_bw()+
  ggtitle("Age Effects Results")+
  ylab("Density")+
  xlab("")+
  scale_fill_manual("Metric",values = renoir_pal[c(1, 6)])+
  scale_color_manual("Metric",values = renoir_pal[c(1, 6)])+
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = "n",
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
age_effect_mse_plot
ggsave(paste0("figures/",modelid,"/age_effect_mse_coverage_density_",modelid,".png"),
      age_effect_mse_plot,
      height = 6,
      width = 10.5)


#####################################################
###
### Period effects bias and coverage
###
#####################################################

period_effect_bias <- vector(mode = "list", length = n_sims)
period_effect_coverage_ind <- vector(mode = "list", length = n_sims)

#removing the couple of seeds that failed for now
#fixing data generating should fix this too
# sims <- 1:n_sims
# sims <- (1:n_sims)[-3]
# sims <- (sims)[-29]


for(i in sims) {
  temp <- fit_sum_ls[[i]][grep("period_effect_surv",rownames(fit_sum_ls[[i]])),]
  period_effect_mn <- temp[, 1]
  period_effect_bias[[i]] <- period_effect_mn - period_effect_true
  coverage_ind <- rep(NA, nT_period)
  for(j in 1:nT_period) {
    coverage_ind[j] <- ifelse(period_effect_true[j]>temp[j,4] & period_effect_true[j] < temp[j,5],1,0)
  }
  period_effect_coverage_ind[[i]] <- coverage_ind
}

period_effect_coverage_mat <- do.call(rbind, period_effect_coverage_ind)
period_indicators <- apply(period_effect_coverage_mat, 2, sum)
period_effect_coverage <- period_indicators / length(sims)

period_effect_bias_df <- do.call(rbind, period_effect_bias)
period_effect_mse <- apply(period_effect_bias_df^2, 2, mean)

plot(period_effect_mse)
plot(period_effect_coverage)


##################################
### Proper plots of period effects
##################################

df_period_effect_out <- data.frame(mse = period_effect_mse,
                                  coverage = period_effect_coverage)

df_period_effect_out$period <- 1:nT_period
df_period_effect_out <- df_period_effect_out %>% pivot_longer(cols = c('mse','coverage')) 
df_period_effect_out$name <- as.factor(df_period_effect_out$name)
levels(df_period_effect_out$name) <- c("Coverage","MSE")

#########################
### period effect mse
#########################
period_mse_coverage_plot <- ggplot(data = df_period_effect_out)+
  geom_point(aes(x = period,y = value,color = name),alpha = .9, size = 1.5)+
  facet_wrap(.~name,scales = "free") +
  theme_bw()+
  ggtitle("Period Effects Results")+
  ylab("")+
  xlab("Period")+
  scale_color_manual(values = renoir_pal[c(1, 6)])+
  # scale_x_continuous(breaks = seq(104,nT_period,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
period_mse_coverage_plot

ggsave(paste0("figures/",modelid,"/period_mse_coverage_plot_",modelid,".png"),
      period_mse_coverage_plot,
      height = 6,
      width = 10.5)

df_period_effect_long <-  data.frame(mse = period_effect_mse,
                                  coverage = period_effect_coverage) %>% pivot_longer(cols = everything())
df_period_effect_long$name <- as.factor(df_period_effect_long$name)
levels(df_period_effect_long$name) <- c("Coverage","MSE")

period_effect_mse_plot <- ggplot(data = df_period_effect_long)+
  geom_density(aes(x = value,color = name,fill = name),alpha = .9)+
  facet_wrap(.~name,scales = "free") +
  theme_bw()+
  ggtitle("Period Effects Results")+
  ylab("Density")+
  xlab("")+
  scale_fill_manual("Metric",values = renoir_pal[c(1, 6)])+
  scale_color_manual("Metric",values = renoir_pal[c(1, 6)])+
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = "n",
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
period_effect_mse_plot

ggsave(paste0("figures/",modelid,"/period_effect_mse_coverage_density_",modelid,".png"),
      period_effect_mse_plot,
      height = 6,
      width = 10.5)


#####################################################
###
### Force of infection effects bias and coverage
###
#####################################################

# foi_age_effect_bias <- vector(mode = "list", length = n_sims)
# foi_age_effect_coverage_ind <- vector(mode = "list", length = n_sims)

#removing the couple of seeds that failed for now
#fixing data generating should fix this too
# sims <- 1:n_sims
# sims <- (1:n_sims)[-3]
# sims <- (sims)[-29]


# for(i in sims) {
#   temp <- fit_sum_ls[[i]][grep("foi_age_effect",rownames(fit_sum_ls[[i]])),]
#   foi_age_effect_mn <- temp[, 1]
#   foi_age_effect_bias[[i]] <- foi_age_effect_mn - foi_age_true
#   coverage_ind <- rep(NA, n_ageclass)
#   for(j in 1:n_ageclass) {
#     coverage_ind[j] <- ifelse(foi_age_true[j]>temp[j,4] & foi_age_true[j] < temp[j,5],1,0)
#   }
#   foi_age_effect_coverage_ind[[i]] <- coverage_ind
# }

# foi_age_effect_coverage_mat <- do.call(rbind, foi_age_effect_coverage_ind)
# foi_age_indicators <- apply(foi_age_effect_coverage_mat, 2, sum)
# foi_age_effect_coverage <- foi_age_indicators / length(sims)

# foi_age_effect_bias_df <- do.call(rbind, foi_age_effect_bias)
# foi_age_effect_mse <- apply(foi_age_effect_bias_df^2, 2, mean)
# foi_age_effect_bias_df <- data.frame(foi_age_effect_bias_df)
# names(foi_age_effect_bias_df) <- c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+")

# foi_age_effect_bias_df_out <- foi_age_effect_bias_df %>% pivot_longer(cols = everything())
# foi_age_effect_bias_df_out$name <- as.factor(foi_age_effect_bias_df_out$name)
# foi_age_effect_bias_df_out$name  <- factor(foi_age_effect_bias_df_out$name, levels = c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+"))

# foi_age_effect_mse_df <- data.frame (mse = as.character(round(foi_age_effect_mse,2)), coverage = as.character(round(foi_age_effect_coverage,2)))
# foi_age_effect_mse_df$ageclass <- c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+")
# foi_age_effect_mse_df$ageclass  <- as.factor(foi_age_effect_mse_df$ageclass)
# foi_age_effect_mse_df$ageclass  <- factor(foi_age_effect_mse_df$ageclass, levels =  c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+"))


# names(foi_age_effect_bias_df_out)
# foi_age_bias_plot <- ggplot(foi_age_effect_bias_df_out,aes(x = name,y=value)) +
#   geom_violin(fill = "lightgrey",color = "lightgrey")+
#   ylim(-1.4,1)+
#   ggtitle("Force of Infection Age Effects Bias, MSE, Coverage")+
#   xlab("Age Class") + ylab("Bias (posterior mean - true value)")+
#   geom_pointrange(stat = "summary",
#                   fun.data = "mean_sdl",
#                   fun.args = list(mult = 2),
#                   color = "darkred",size = 1.5)+
#   geom_hline(yintercept =  0, linetype = "dashed", size = 2) + 
#   geom_text(data = foi_age_effect_mse_df, 
#              aes(x = ageclass,y = .98,
#              label = paste0("MSE:",mse)))+
#   geom_text(data = foi_age_effect_mse_df, 
#              aes(x = ageclass,y = .9,
#              label = paste0("Cover:",coverage)))+
#   theme_bw()

# foi_age_bias_plot

# ggsave(paste0("figures/",modelid,"/foi_age_bias_plot_",modelid,".png"),foi_age_bias_plot,
#         height = 6,width = 10.5)

# plot(foi_age_effect_mse)
# plot(foi_age_effect_coverage)

#####################################################
###
### Outputting aggregated results
###
#####################################################


age_effect = c(sum(age_effect_mse),mean(age_effect_coverage))
period_effect = c(sum(period_effect_mse), mean(period_effect_coverage))


survival_params_tab <- data.frame(beta0_tab,age_effect = c(sum(age_effect_mse),mean(age_effect_coverage)),
period_effect = c(sum(period_effect_mse), mean(period_effect_coverage)))
survival_params_tab <- round(survival_params_tab,3)

write.csv(survival_params_tab, file =paste0("results/",modelid,"/survival_params_tab_",modelid,".csv"))
write.csv(survival_params_tab, file =paste0("figures/",modelid,"/survival_params_tab_",modelid,".csv"))


#####################################################
###
### Survival Probabilities
###
#####################################################



#####################################################
###
### sn_sus mse, coverage, and bias plots
###
#####################################################

sn_sus_bias <- vector(mode = "list", length = n_sims)
sn_sus_coverage_ind <- vector(mode = "list", length = n_sims)

for(i in sims) {
  temp <- fit_sum_ls[[i]][grep("sn_sus",rownames(fit_sum_ls[[i]])),]
  sn_sus_mn <- temp[, 1]
  sn_sus_bias[[i]] <- sn_sus_mn - sn_sus_true
  coverage_ind <- rep(NA, n_ageclass)
  for(j in 1:length(sn_sus_true)) {
    coverage_ind[j] <- ifelse(sn_sus_true[j]>temp[j,4] & sn_sus_true[j] < temp[j,5],1,0)
  }
  sn_sus_coverage_ind[[i]] <- coverage_ind
}

sn_sus_coverage_mat <- do.call(rbind, sn_sus_coverage_ind)
sn_sus_indicators <- apply(sn_sus_coverage_mat, 2, sum)
sn_sus_coverage <- sn_sus_indicators / length(sims)

sn_sus_bias_df <- do.call(rbind, sn_sus_bias)
sn_sus_mse <- apply(sn_sus_bias_df^2, 2, mean)

sn_sus_df_out <- data.frame(mse = sn_sus_mse,
                            coverage = sn_sus_coverage,
                            year = rep(1:5,each = 7),
                            age = rep(1:7,5))
sn_sus_df_out$age <- factor(sn_sus_df_out$age)
levels(sn_sus_df_out$age) <- c("Fawn","1.5","2.5","3.5","4.5","5.5","6.5+")


sn_sus_mse_plot <- ggplot(sn_sus_df_out,aes(x = year))+
  geom_point(aes(y = mse,color = age),size = 3)+
  facet_wrap(~age,nrow = 2)+
  theme_bw()+
  ggtitle("Susceptible Annual Survival MSE Results")+
  ylab("MSE")+
  xlab("Year")+
  scale_color_manual("Age",values=renoir_pal[1:7])+
  # scale_color_manual(values = renoir_pal[c(1, 6)])+
  # scale_x_continuous(breaks = seq(104,nT_period,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )

ggsave(paste0("figures/",modelid,"/sn_sus_mse_plot_",modelid,".png"),
        sn_sus_mse_plot,
        height = 6,width = 10.5)



sn_sus_cover_plot <- ggplot(sn_sus_df_out,aes(x = year))+
  geom_point(aes(y = coverage,color = age),size = 3)+
  facet_wrap(~age,nrow = 2)+
  theme_bw()+
  ggtitle("Susceptible Annual Survival Coverage Results")+
  ylab("Coverage")+
  xlab("Year")+
  ylim(0,1)+
  scale_color_manual("Age",values=renoir_pal[1:7])+
  # scale_color_manual(values = renoir_pal[c(1, 6)])+
  # scale_x_continuous(breaks = seq(104,nT_period,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )

ggsave(paste0("figures/",modelid,"/sn_sus_cover_plot_",modelid,".png"),
        sn_sus_cover_plot,
        height = 6,width = 10.5)

###
### bias
###

sn_sus_bias_df <- data.frame(sn_sus_bias_df) %>% pivot_longer(cols = everything())
sn_sus_bias_df$age <- as.factor(rep(rep(1:7,5),length(sims)))
levels(sn_sus_bias_df$age) <- c("Fawn","1.5","2.5","3.5","4.5","5.5","6.5+")
sn_sus_bias_df$year <- as.factor(rep(rep(1:5,each = 7),length(sims)))

sn_sus_bias_plot <- ggplot(sn_sus_bias_df,aes(x = year,y=value,color = age,fill = age)) +
  geom_violin(alpha = .6)+
  facet_wrap(~age,nrow = 2)+
  ggtitle("Susceptible Annual Survival Probability Bias")+
  xlab("Year") + ylab("Bias (posterior mean - true value)")+
  geom_pointrange(stat = "summary",
                  fun.data = "mean_sdl",
                  fun.args = list(mult = 2),size = 1.5)+
  geom_hline(yintercept =  0, linetype = "dotted", size = 1) + 
    scale_fill_manual("Age",values=renoir_pal[1:7])+
    scale_color_manual("Age",values=renoir_pal[1:7])+
  # geom_text(data = sn_sus_mse_df, 
  #            aes(x = ageclass,y = .98,
  #            label = paste0("MSE:",mse)))+
  # geom_text(data = sn_sus_mse_df, 
  #            aes(x = ageclass,y = .9,
  #            label = paste0("Cover:",coverage)))+
  theme_bw()+theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
sn_sus_bias_plot
ggsave(paste0("figures/",modelid,"/sn_sus_bias_plot_",modelid,".png"),sn_sus_bias_plot,
        height = 6,width = 10.5)

#####################################################
###
### sn_inf mse, coverage, and bias plots
###
#####################################################

sn_inf_bias <- vector(mode = "list", length = n_sims)
sn_inf_coverage_ind <- vector(mode = "list", length = n_sims)

#removing the couple of seeds that failed for now
#fixing data generating should fix this too
# sims <- 1:n_sims
# sims <- (1:n_sims)[-3]
# sims <- (sims)[-29]

for(i in sims) {
  temp <- fit_sum_ls[[i]][grep("sn_inf",rownames(fit_sum_ls[[i]])),]
  sn_inf_mn <- temp[, 1]
  sn_inf_bias[[i]] <- sn_inf_mn - sn_inf_true
  coverage_ind <- rep(NA, n_ageclass)
  for(j in 1:length(sn_inf_true)) {
    coverage_ind[j] <- ifelse(sn_inf_true[j]>temp[j,4] & sn_inf_true[j] < temp[j,5],1,0)
  }
  sn_inf_coverage_ind[[i]] <- coverage_ind
}

sn_inf_coverage_mat <- do.call(rbind, sn_inf_coverage_ind)
sn_inf_indicators <- apply(sn_inf_coverage_mat, 2, sum)
sn_inf_coverage <- sn_inf_indicators / length(sims)

sn_inf_bias_df <- do.call(rbind, sn_inf_bias)
sn_inf_mse <- apply(sn_inf_bias_df^2, 2, mean)

sn_inf_df_out <- data.frame(mse = sn_inf_mse,
                            coverage = sn_inf_coverage,
                            year = rep(1:5,each = 7),
                            age = rep(1:7,5))
sn_inf_df_out$age <- factor(sn_inf_df_out$age)
levels(sn_inf_df_out$age) <- c("Fawn","1.5","2.5","3.5","4.5","5.5","6.5+")


sn_inf_mse_plot <- ggplot(sn_inf_df_out,aes(x = year))+
  geom_point(aes(y = mse,color = age),size = 3)+
  facet_wrap(~age,nrow = 2)+
  theme_bw()+
  ggtitle("Infected Annual Survival MSE Results")+
  ylab("MSE")+
  xlab("Year")+
  scale_color_manual("Age",values=renoir_pal[1:7])+
  # scale_color_manual(values = renoir_pal[c(1, 6)])+
  # scale_x_continuous(breaks = seq(104,nT_period,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )

ggsave(paste0("figures/",modelid,"/sn_inf_mse_plot_",modelid,".png"),
        sn_inf_mse_plot,
        height = 6,width = 10.5)



sn_inf_cover_plot <- ggplot(sn_inf_df_out,aes(x = year))+
  geom_point(aes(y = coverage,color = age),size = 3)+
  facet_wrap(~age,nrow = 2)+
  theme_bw()+
  ggtitle("Infected Annual Survival Coverage Results")+
  ylab("Coverage")+
  xlab("Year")+
  ylim(0,1)+
  scale_color_manual("Age",values=renoir_pal[1:7])+
  # scale_color_manual(values = renoir_pal[c(1, 6)])+
  # scale_x_continuous(breaks = seq(104,nT_period,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )

ggsave(paste0("figures/",modelid,"/sn_inf_cover_plot_",modelid,".png"),
        sn_inf_cover_plot,
        height = 6,width = 10.5)

###
### bias
###

sn_inf_bias_df <- data.frame(sn_inf_bias_df) %>% pivot_longer(cols = everything())
sn_inf_bias_df$age <- as.factor(rep(rep(1:7,5),length(sims)))
levels(sn_inf_bias_df$age) <- c("Fawn","1.5","2.5","3.5","4.5","5.5","6.5+")
sn_inf_bias_df$year <- as.factor(rep(rep(1:5,each = 7),length(sims)))
sn_inf_bias_df <- sn_inf_bias_df %>% filter(age != "Fawn")

sn_inf_bias_plot <- ggplot(sn_inf_bias_df,aes(x = year,y=value,color = age,fill = age)) +
  geom_violin(alpha = .6)+
  facet_wrap(~age,nrow = 2)+
  ggtitle("Infected Annual Survival Probability Bias")+
  xlab("Year") + ylab("Bias (posterior mean - true value)")+
  geom_pointrange(stat = "summary",
                  fun.data = "mean_sdl",
                  fun.args = list(mult = 2),size = 1.5)+
  geom_hline(yintercept =  0, linetype = "dotted", size = 1) + 
    scale_fill_manual("Age",values=renoir_pal[1:7])+
    scale_color_manual("Age",values=renoir_pal[1:7])+
  # geom_text(data = sn_inf_mse_df, 
  #            aes(x = ageclass,y = .98,
  #            label = paste0("MSE:",mse)))+
  # geom_text(data = sn_inf_mse_df, 
  #            aes(x = ageclass,y = .9,
  #            label = paste0("Cover:",coverage)))+
  theme_bw()+theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    strip.text = element_text(size = 16),
                    legend.position = 'n',
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14)
              )
sn_inf_bias_plot
ggsave(paste0("figures/",modelid,"/sn_inf_bias_plot_",modelid,".png"),sn_inf_bias_plot,
        height = 6,width = 10.5)


#####################################################
###
### Outputting aggregated results
###
#####################################################


sink(paste0("results/",modelid,"/results_simulation_",modelid,".txt"))
cat(paste0("beta0_survival_sus_mse","\n", beta0_survival_sus_mse,"\n"))
cat(paste0("beta0_survival_suscoverage","\n", beta0_survival_sus_coverage,"\n"))
cat(paste0("mean runtime: ",runtime_mn,"\n"))
cat("age_effect_mse: \n", age_effect_mse,"\n")
cat("age_effect_coverage: \n", age_effect_coverage,"\n")
cat("period_effect_mse: \n", period_effect_mse,"\n")
cat("period_effect_coverage: \n", period_effect_coverage,"\n")
cat("effective sample size: \n")
# print(ess_ls)
sink()

save(beta0_survival_sus_bias,file=paste0("results/",modelid,"/beta0_survival_sus_bias_",modelid,".Rdata"))
save(beta0_survival_sus_mse,file=paste0("results/",modelid,"/beta0_survival_sus_mse_",modelid,".Rdata"))
save(beta0_survival_sus_coverage,file=paste0("results/",modelid,"/beta0_survival_sus_coverage_",modelid,".Rdata"))
save(runtime_mn,file=paste0("results/",modelid,"/runtime_mn_",modelid,".Rdata"))
save(age_effect_bias,file=paste0("results/",modelid,"/age_effect_bias_",modelid,".Rdata"))
save(age_effect_mse,file=paste0("results/",modelid,"/age_effect_mse_",modelid,".Rdata"))
save(age_effect_coverage,file=paste0("results/",modelid,"/age_effect_coverage_",modelid,".Rdata"))
save(period_effect_bias,file=paste0("results/",modelid,"/period_effect_bias_",modelid,".Rdata"))
save(period_effect_mse,file=paste0("results/",modelid,"/period_effect_mse_",modelid,".Rdata"))
save(period_effect_coverage,file=paste0("results/",modelid,"/period_effect_coverage_",modelid,".Rdata"))
save(ess_ls,file=paste0("results/",modelid,"/ess_ls_",modelid,".Rdata"))