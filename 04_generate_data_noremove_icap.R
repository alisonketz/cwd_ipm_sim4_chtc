##################################################################
###
### run the function to generate data
###
##################################################################

source("03_fun_generate_data.R")

dat <- ageperiod_surv_foi_sim_data(
            beta0_survival_sus = beta0_survival_sus,
            beta0_survival_inf = beta0_survival_inf,
            foi_age_effect = foi_age_effect,
            age_effect = age_effect_true,
            period_effect = period_effect_true,
            nT_age = nT_age,
            nT_period = nT_period,
            processnum = processnum
            )

n_fit <- nrow(dat$df_fit)
n_fit_sus_cens_posttest <- nrow(dat$d_fit_sus_cens_posttest)
n_fit_sus_cens_postno <- nrow(dat$d_fit_sus_cens_postno)
n_fit_sus_mort_posttest <- nrow(dat$d_fit_sus_mort_posttest)
n_fit_sus_mort_postno <- nrow(dat$d_fit_sus_mort_postno)
n_fit_icap_cens <- nrow(dat$d_fit_icap_cens)
n_fit_icap_mort <- nrow(dat$d_fit_icap_mort)
n_fit_rec_neg_cens_posttest <- nrow(dat$d_fit_rec_neg_cens_posttest)
n_fit_rec_neg_cens_postno <- nrow(dat$d_fit_rec_neg_cens_postno)
n_fit_rec_neg_mort <- nrow(dat$d_fit_rec_neg_mort)
n_fit_rec_pos_cens <- nrow(dat$d_fit_rec_pos_cens)
n_fit_rec_pos_mort <- nrow(dat$d_fit_rec_pos_mort)
n_fit_idead <- nrow(dat$d_fit_idead)

### test that these data combined are the same
### dimension as the overall generated data
# n_fit_sus_cens_posttest+
# n_fit_sus_cens_postno+
# n_fit_sus_mort_posttest+
# n_fit_sus_mort_postno+
# n_fit_icap_cens+
# n_fit_icap_mort+
# n_fit_rec_neg_cens_posttest+
# n_fit_rec_neg_cens_postno+
# n_fit_rec_neg_mort+
# n_fit_rec_pos_cens+
# n_fit_rec_pos_mort+
# n_fit_idead

datatypes_out <- data.frame(datatype = c("overall",
"sus_cens_posttest",
"sus_cens_postno",
"sus_mort_posttest",
"sus_mort_postno",
"icap_cens",
"icap_mort",
"rec_neg_cens_posttest",
"rec_neg_cens_postno",
"rec_neg_mort",
"rec_pos_cens",
"rec_pos_mort",
"idead"),
number_samples = c(n_fit,
n_fit_sus_cens_posttest,
n_fit_sus_cens_postno,
n_fit_sus_mort_posttest,
n_fit_sus_mort_postno,
n_fit_icap_cens,
n_fit_icap_mort,
n_fit_rec_neg_cens_posttest,
n_fit_rec_neg_cens_postno,
n_fit_rec_neg_mort,
n_fit_rec_pos_cens,
n_fit_rec_pos_mort,
n_fit_idead))

write.csv(datatypes_out, file = paste0("datatypes_out_",processnum,".csv"))


##############################################################
### there are too many deer infected at to capture
### because we aren't accounting for 
### disease-associated mortality in the generating function
##############################################################
# icap_eval_df <- data.frame(left_age = dat$left_age[dat$pos1==TRUE])
# icap_eval_plot <- ggplot(data = icap_eval_df) +
#       geom_histogram(aes(x = left_age),bins = 50) +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )


# ggsave(paste0("figures/icap_eval_left_age.png"),
#       icap_eval_plot,
#       height = 6,
#       width = 8)

# icap_eval_df$weights <- icap_eval_df$left_age/sum(icap_eval_df$left_age)
# icap_eval_df$weights2 <- (1 - 1/icap_eval_df$left_age)/sum((1 - 1/icap_eval_df$left_age))

# icap_weights_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# icap_weights2_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights2), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# ggsave(paste0("figures/icap_eval_weights.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)



# ggsave(paste0("figures/icap_eval_weights2.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)


# n_remove <- sum(dat$pos1) - round(sum(dat$pos1)/5)
# indx_remove <- sample(which(dat$pos1==TRUE),n_remove,replace=FALSE) 
# n_ind <- dat$n_ind - n_remove
# n_ind <- dat$n_ind



#####################################################
###
### reformatting to fit in survival only type model
###
#####################################################

d_surv <- dat$df_fit
n_surv <- nrow(d_surv)

d_surv$cwd_mort[d_surv$cwd_cap == 1] <- 1
sum(d_surv$cwd_mort)/1250

fast <- d_surv$id[d_surv$r_age == d_surv$e_age]
d_fit <- data.frame(matrix(NA,nr = n_surv + sum(d_surv$rt_censor == 0), ncol = ncol(d_surv)))
names(d_fit) <- names(d_surv)
d_fit[1:n_surv,] <- d_surv
n_fit <- nrow(d_fit)
d_temp <- d_surv[d_surv$rt_censor == 0,]
d_temp$rt_censor  <- 1
d_fit[(n_surv + 1):n_fit,] <- d_temp
names(d_fit) <- names(d_surv)
d_fit <- d_fit[order(d_fit$id),]
d_fit$left_age <- d_fit$e_age
d_fit$left_age[d_fit$rt_censor == 0] <- d_fit$r_age[d_fit$rt_censor == 0]
d_fit$left_period <- d_fit$e_period
d_fit$left_period[d_fit$rt_censor == 0] <- d_fit$r_period[d_fit$rt_censor == 0]

d_fit$right_age <- d_fit$r_age
d_fit$right_age[d_fit$rt_censor == 0] <- d_fit$s_age[d_fit$rt_censor == 0]
d_fit$right_period <- d_fit$r_period
d_fit$right_period[d_fit$rt_censor == 0] <- d_fit$s_period[d_fit$rt_censor == 0]

d_fit$age2date <- d_fit$left_period - d_fit$left_age
rm <- which(d_fit$id %in% fast & d_fit$rt_censor == 1)
d_fit <- d_fit[-rm,]
n_fit <- nrow(d_fit)

low_pos_mort <- d_surv$id[which(d_surv$cwd_mort == 1)]
d_fit_pos <- d_fit[d_fit$id %in% low_pos_mort,]
n_fit_pos <- nrow(d_fit_pos)
d_fit_neg <- d_fit[!(d_fit$id %in% low_pos_mort),]
n_fit_neg <- nrow(d_fit_neg)


sum(1-d_fit$rt_censor)
d_fit_pos <- d_fit[d_fit$cwd_mort == 1,]
head(d_fit[d_fit$cwd_cap ==1,])
d_fit_neg <- d_fit[d_fit$cwd_mort == 0,]