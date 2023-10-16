
##################################################################
###
###  Function to calculate the
###  probability in the likelihood
###
##################################################################

dSurv <- nimbleFunction(
  run = function(x= integer(0),
                 n_fit = integer(0),
                 censored = double(1),
                 left = double(1),
                 right = double(1),
                 age_effect = double(1),
                 period_effect = double(1),
                 age2date = double(1),
                 nT_age = double(0),
                 beta0 = double(0),
                 log = double(0)
  ){
    sumllik <- 0
    for (j in 1:n_fit) {
      UCH <- nimNumeric(nT_age)
      for (k in left[j]:(right[j] - 1)) {
        UCH[k] <- exp(beta0 +
                      age_effect[k] +
                      period_effect[k + age2date[j]]
                      )
      }
      p <- exp(-sum(UCH[left[j]:(right[j] - 1)]))
      sumllik <- sumllik + dbinom(censored[j], 1, p, log=TRUE)
    if(is.na(sumllik)){stop("ack")}
    }
    returnType(double(0))
    if (log) {
        return(sumllik)
    } else {
        return(exp(sumllik))
    }
})


nimble::registerDistributions(list(
    dSurv = list(
        BUGSdist = "dSurv(n_fit,censored,left,right,age2date,age_effect,period_effect,nT_age,beta0)",
        types = c(
            "value= integer(0)",
            "n_fit = integer(0)",
            "censored = double(1)",
            "left = double(1)",
            "right = double(1)",
            "age2date = double(1)",
            "age_effect = double(1)",
            "period_effect = double(1)",
            "nT_age = double(0)",
            "beta0 = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign("dSurv",
    dSurv,
    envir = .GlobalEnv
)







# #######################################################################
# ###
# ###   Likelihoods for each kind of data listed below
# ###
# #######################################################################
# # d_fit_icap_mort
# # d_fit_icap_cens
# # d_fit_idead
# # d_fit_sus_mort_posttest
# # d_fit_sus_cens_postno
# # d_fit_sus_cens_posttest
# # d_fit_sus_mort_postno
# # d_fit_rec_neg_cens_posttest
# # d_fit_rec_neg_cens_postno
# # d_fit_rec_neg_mort
# # d_fit_rec_pos_cens
# # d_fit_rec_pos_mort

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer mortalities for radio marked deer that
# ###   enter the study as test positive at capture
# ###
# ###   d_fit_icap_mort
# ###
# ###   Overleaf Equation (12)
# ###
# #######################################################################

# dIcapMort <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    age2date = double(1),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    log = double()) {

#         # start the loop through individuals
#         sumllik <- 0
#         for (i in 1:n_samples) {
#             lam_inf <- 0
#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################

#             if(r[i]>e[i]) {
#                     for (j in e[i]:(r[i] - 1)) {
#                     lam_inf <- lam_inf +
#                                exp(beta0_inf +
#                                    age_effect_surv[j] +
#                                    period_effect_surv[age2date[i] + j]
#                                    )
#                     }
#                 }
#             lam_inf_s <- sum(exp(beta0_inf +
#                 age_effect_surv[r[i]:(s[i] - 1)] +
#                 period_effect_surv[age2date[i] + r[i]:(s[i] - 1)]))
                
#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             sumllik <- sumllik - lam_inf +
#                        log(1 - exp(-lam_inf_s))
#         }#end n_samples

#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
# })


# nimble::registerDistributions(list(
#     dIcapMort = list(
#         BUGSdist = "dIcapMort(n_samples,e,r,s,age2date,beta0_inf,age_effect_surv,period_effect_surv)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "age2date = double(1)",
#             "beta0_inf = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# assign("dIcapMort", dIcapMort, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dIcapMort(
# #         x = 1,
# #         n_samples = nrow(d_fit_icap_mort),
# #         e = d_fit_icap_mort$e_age,
# #         r = d_fit_icap_mort$r_age,
# #         s = d_fit_icap_mort$s_age,
# #         age2date = d_fit_icap_mort$age2date,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_survival_test,
# #         period_effect_surv = period_effect_surv,
# #         log = TRUE
# #         )
# # (end <- Sys.time() - starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer mortalities for radio marked deer that
# ###   enter the study as test positive at capture
# ###
# ###   d_fit_icap_cens
# ###
# ###   Overleaf Equation (10)
# ###
# #######################################################################

# dIcapCens <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    age2date = double(1),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_inf <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################

#             # recalculate lam_inf to sum over e to r-1
#             for (j in e[i]:(r[i] - 1)) {
#                 lam_inf <- lam_inf +
#                             exp(beta0_inf +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j])
#             }

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_inf
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dIcapCens = list(
#         BUGSdist = "dIcapCens(n_samples,e,r,age2date,beta0_inf,age_effect_surv,period_effect_surv)",
#         types = c(
#             "value=integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "age2date = double(1)",
#             "beta0_inf = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# # for a user-defined distribution
# assign("dIcapCens", dIcapCens, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dIcapCens(
# #         x = 1,
# #         n_samples = nrow(d_fit_icap_cens),
# #         e = d_fit_icap_cens$e_age,
# #         r = d_fit_icap_cens$r_age,
# #         age2date = icap_cens_age2date,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_survival_test,
# #         period_effect_surv = period_effect_surv,
# #         log = TRUE
# #         )
# # (endtime <- Sys.time() - starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer mortalities for radio marked deer that
# ###   enter the study as test negative at capture
# ###
# ###   d_fit_idead
# ###
# ###   Overleaf Equation (24)
# ###
# #######################################################################

# dNegCapPosMort <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {

#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
#             lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             # survival hazard for susceptible deer
#             lam_sus[e[i]:(s[i] - 1)] <- exp(beta0_sus +
#                 age_effect_surv[e[i]:(s[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (s[i] - 1 + age2date[i])])
#             # survival hazard while infected
#             lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
#                 age_effect_surv[e[i]:(s[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (s[i] - 1 + age2date[i])])
#             # force of infection infection hazard
#             lam_foi[e[i]:(s[i] - 1)] <- exp(foi_age_effect[age_lookup[e[i]:(s[i] - 1)]])
        

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             lik_temp <- lik_temp +
#                         (1 - exp(-lam_foi[e[i]])) *
#                         exp(-sum(lam_inf[e[i]:(r[i] - 1)]))
#             if((r[i] - e[i])>2) {
#                 for (k in (e[i] + 1):(r[i] - 1)) {
#                         lik_temp <- lik_temp +
#                                     (1 - exp(-lam_foi[k])) *
#                                     exp(-sum(lam_inf[(k):(r[i] - 1)])) *
#                                     exp(-sum(lam_sus[e[i]:(k - 1)])) *
#                                     exp(-sum(lam_foi[e[i]:(k - 1)]))
#                 }
#             }
#             lik_temp <- lik_temp +
#                         (1 - exp(-lam_foi[s[i] - 1])) *
#                         exp(-sum(lam_sus[e[i]:((s[i] - 1) - 1)])) *
#                         exp(-sum(lam_foi[e[i]:((s[i] - 1) - 1)]))

#             sumllik <- sumllik +
#                        log(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) +
#                        log(lik_temp)
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dNegCapPosMort = list(
#         BUGSdist = "dNegCapPosMort(n_samples,e,r,s,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "beta0_inf = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# assign("dNegCapPosMort", dNegCapPosMort, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dNegCapPosMort(
# #         x = 1,
# #         n_samples = nrow(df_fit_idead),
# #         e = df_fit_idead$e_age,
# #         r = df_fit_idead$r_age,
# #         s = df_fit_idead$s_age,
# #         age2date = df_fit_idead$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = age_foi,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected radio-marked deer mortalities:
# ###   test neg at cap and tested mort
# ###
# ###   d_fit_sus_mort_posttest
# ###
# ###   Overleaf Equation (6)
# ###
# #######################################################################

# dSusMortTest <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = double(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0
#             lam_sus <- 0
#             lam_sus_s <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             for (j in r[i]:(s[i] - 1)) {
#                 # sum up surv hazard from r to s-1
#                 lam_sus_s <- lam_sus_s +
#                                 exp(beta0_sus +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j])
#             }

#             for (j in e[i]:(s[i] - 1)) {
#                 #  up foi hazard from 1  to j
#                 lam_foi <- lam_foi +
#                             exp(foi_age_effect[age_lookup[j]])
#                 if (j < r[i]) {
#                 lam_sus <- lam_sus +
#                             exp(beta0_sus +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j])
#                 }
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             sumllik <- sumllik +
#                        log(1 - exp(-lam_sus_s)) -
#                        lam_sus -
#                        lam_foi
#         }

#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dSusMortTest = list(
#         BUGSdist = "dSusMortTest(n_samples,e,r,s,age2date,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = double(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# # ###for a user-defined distribution
# assign("dSusMortTest", dSusMortTest, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dSusMortTest(
# #         x = 1,
# #         n_samples = nrow(d_fit_sus_mort_posttest),
# #         e = d_fit_sus_mort_posttest$e_age,
# #         r = d_fit_sus_mort_posttest$r_age,
# #         s = d_fit_sus_mort_posttest$s_age,
# #         age2date = sus_mort_posttest_age2date,
# #         beta0_sus = beta0_survival_sus,
# #         age_effect_surv = age_effect_survival_test,
# #         period_effect_surv = period_effect_surv,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   Uninfected radio-marked deer right censored:
# ###   Test neg at cap and censoring
# ###
# ###   d_fit_sus_cens_postno
# ###
# ###   Overleaf Equation (4)
# ###
# #######################################################################

# dSusCensNo <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(r[i] - 1)
#             lam_sus <- nimNumeric(r[i] - 1)
#             lam_inf <- nimNumeric(r[i]- 1)
#             lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             lam_foi[e[i]:(r[i] - 1)] <- exp(foi_age_effect[age_lookup[e[i]:(r[i]-1)]])

#             lam_sus[e[i]:(r[i] - 1)] <- exp(beta0_sus +
#                                         age_effect_surv[e[i]:(r[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (r[i] - 1 + age2date[i])])
#             lam_inf[e[i]:(r[i] - 1)] <- exp(beta0_inf +
#                     age_effect_surv[e[i]:(r[i] - 1)] +
#                     period_effect_surv[(age2date[i] + e[i]):
#                                         (age2date[i] + r[i] - 1)])

#             lik_temp <- (1 - exp(-lam_foi[e[i]])) *
#                         exp(-sum(lam_inf[e[i]:(r[i] - 1)]))

#             if ((r[i] - e[i]) > 1) {
#                 for(k in (e[i] + 1):(r[i] - 1)) {
#                     lik_temp <- lik_temp +
#                                 (1 - exp(-lam_foi[k])) *
#                                 exp(-sum(lam_inf[k:(r[i] - 1)])) *
#                                 exp(-sum(lam_foi[e[i]:(k - 1)])) *
#                                 exp(-sum(lam_sus[e[i]:(k - 1)])) 
#                 }
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik +
#                        log(
#                         exp(-sum(lam_foi[e[i]:(r[i] - 1)])) *
#                         exp(-sum(lam_sus[e[i]:(r[i] - 1)])) +
#                         lik_temp
#                        )
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dSusCensNo = list(
#         BUGSdist = "dSusCensNo(n_samples,e,r,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value=integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "beta0_inf = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# assign("dSusCensNo", dSusCensNo, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dSusCensNo(
# #         x = 1,
#         # n_samples = nrow(d_fit_sus_cens_postno),
#         # e = d_fit_sus_cens_postno$e_age,    
#         # r = d_fit_sus_cens_postno$r_age,   
#         # age2date = sus_cens_postno_age2date, 
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_survival_test,
# #         period_effect_surv = period_effect_surv,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   Uninfected radio-marked deer right censor:
# ###   Test neg at cap and censoring
# ###
# ###   d_fit_sus_cens_posttest
# ###   Overleaf Equation 2
# ###
# #######################################################################

# dSusCensTest <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0),
#                    e = double(1), 
#                    r = double(1), 
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {

#         # starttime the loop through individuals
#         sumllik <- 0
#         for (i in 1:n_samples) {
#             # intitialize vectors
#             lam_foi <- 0
#             lam_sus <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################

#             for (j in e[i]:(r[i] - 1)) {
#                 # sum up foi
#                 lam_foi <- lam_foi +
#                     exp(foi_age_effect[age_lookup[j]])

#                 lam_sus <- lam_sus +
#                         exp(beta0_sus +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j])
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             sumllik <- sumllik - lam_sus - lam_foi

#         } # end loop over individuals

#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dSusCensTest = list(
#         BUGSdist = "dSusCensTest(n_samples,e,r,age2date,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))


# ### for a user-defined distribution
# assign("dSusCensTest", dSusCensTest, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dSusCensTest(
# #         x = 1,
# #         n_samples = n_fit_sus_cens_posttest,
# #         e = dat$d_fit_sus_cens_posttest$e_age,
# #         r = dat$d_fit_sus_cens_posttest$r_age,
# #         age2date = dat$d_fit_sus_cens_posttest$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = age_foi,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test



# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected radio-marked deer mortalities:
# ###   test neg at cap and no test at mortality, no recap
# ###
# ###   d_fit_sus_mort_postno
# ###
# ###   Overleaf Equation (8)
# ###
# #######################################################################

# dSusMortNoTest <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double(0)) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
#             lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             ### survival hazard for susceptible deer
#             lam_sus[e[i]:(s[i] - 1)] <- exp(beta0_sus +
#                 age_effect_surv[e[i]:(s[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (s[i] - 1 + age2date[i])])

#             ### survival hazard while infected
#             lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
#                 age_effect_surv[e[i]:(s[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (s[i] - 1 + age2date[i])])

#             ### force of infection infection hazard
#             lam_foi[e[i]:(s[i] - 1)] <- 
#                 exp(foi_age_effect[age_lookup[e[i]:(s[i] - 1)]])

#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             lik_temp <- lik_temp +
#                         (1 - exp(-lam_foi[e[i]])) *
#                         (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) *
#                         exp(-sum(lam_inf[e[i]:(r[i] - 1)]))
#             if((r[i] - e[i]) > 2) {
#                 for(k in (e[i] + 1):(r[i] - 1)) {
#                     lik_temp <- lik_temp + 
#                                 (1 - exp(-lam_foi[k])) *
#                                 (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) *
#                                 exp(-sum(lam_foi[e[i]:(k - 1)])) *
#                                 exp(-sum(lam_sus[e[i]:(k - 1)])) *
#                                 exp(-sum(lam_inf[k:(r[i] - 1)]))
#                 }
#                 lik_temp <- lik_temp +
#                             (1 - exp(-lam_foi[s[i]-1])) * 
#                             (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) *
#                             exp(-sum(lam_foi[e[i]:(r[i] - 1)])) *
#                             exp(-sum(lam_sus[e[i]:(r[i] - 1)]))
#                 sumllik <- sumllik +
#                        log(
#                         exp(-sum(lam_foi[e[i]:(s[i] - 1)])) *
#                         exp(-sum(lam_sus[e[i]:(r[i] - 1)])) *
#                         (1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))) +
#                         lik_temp
#                        )#endloglik
#             } else {
#                 lik_temp <- lik_temp + 
#                                 (1 - exp(-lam_foi[e[i]])) *
#                                 (1 - exp(-sum(lam_inf[e[i]:(s[i] - 1)])))

#                 sumllik <- sumllik +
#                        log(
#                         exp(-sum(lam_foi[e[i]:(s[i] - 1)])) *
#                         (1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)])))) +
#                         lik_temp

#             }
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         }
#     }
# )

# nimble::registerDistributions(list(
#     dSusMortNoTest = list(
#         BUGSdist = "dSusMortNoTest(n_samples,e,r,s,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "age2date = double(1)",
#             "beta0_inf = double(0)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# assign("dSusMortNoTest", dSusMortNoTest, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dSusMortNoTest(
# #         x = 1,
# #         n_samples = n_fit_sus_mort_postno,
# #         e = dat$d_fit_sus_mort_postno$e_age,
# #         r = dat$d_fit_sus_mort_postno$r_age,
# #         s = dat$d_fit_sus_mort_postno$s_age,
# #         age2date = dat$d_fit_sus_mort_postno$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected deer that were test neg at capture,
# ###   then test negative at recap, that are right censored,
# ###   and have been tested post censoring
# ###
# ###   d_fit_rec_neg_cens_posttest
# ###
# ###   Overleaf Equation (14)
# ###
# #######################################################################

# dRecNegCensTest <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0
#             lam_sus <- 0
#             # lam_inf <- 0
#             # lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             for (j in e[i]:(r[i] - 1)) {
#                 ### survival hazard for susceptible deer
#                 lam_sus <- lam_sus + exp(beta0_sus +
#                     age_effect_surv[j] +
#                     period_effect_surv[j + age2date[i]])
#                 ### force of infection infection hazard
#                 lam_foi <- lam_foi + exp(foi_age_effect[age_lookup[j]])
#             }

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_sus - lam_foi
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecNegCensTest = list(
#         BUGSdist = "dRecNegCensTest(n_samples,e,r,age2date,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ## for a user-defined distribution
# assign("dRecNegCensTest", dRecNegCensTest, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dRecNegCensTest(
# #         x = 1,
# #         n_samples = n_fit_rec_neg_cens_posttest,
# #         e = dat$d_fit_rec_neg_cens_posttest$e_age,
# #         r = dat$d_fit_rec_neg_cens_posttest$r_age,
# #         age2date = dat$d_fit_rec_neg_cens_posttest$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected deer that were test neg at capture,
# ###   then test negative at recap, that are right censored,
# ###   and have not been tested post censoring
# ###
# ###   d_fit_rec_neg_cens_postno
# ###
# ###   Overleaf Equation (16)
# ###
# #######################################################################

# dRecNegCensPostNo <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    dn1 = double(1), # right of last test negative
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(r[i], init = FALSE)
#             lam_sus <- nimNumeric(r[i], init = FALSE)
#             lam_inf <- nimNumeric(r[i], init = FALSE)
#             # lik_temp <- 0

#             #############################################
#             ### preliminary hazards for the likelihood
#             #############################################
#             ### survival hazard for susceptible deer
#             lam_sus[e[i]:(r[i] - 1)] <- exp(beta0_sus +
#                 age_effect_surv[e[i]:(r[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (r[i] - 1 + age2date[i])])
#             ### survival hazard while infected
#             lam_inf[(dn1[i]):(r[i] - 1)] <- exp(beta0_inf +
#                 age_effect_surv[(dn1[i]):(r[i] - 1)] +
#                 period_effect_surv[(dn1[i] + age2date[i]):
#                                     (r[i] - 1 + age2date[i])])
#             ### force of infection infection hazard
#             lam_foi[e[i]:(r[i] - 1)] <-
#                 exp(foi_age_effect[age_lookup[e[i]:(r[i] - 1)]])

#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             if((r[i] - dn1[i]) > 2) {
#                 lik_temp <- (1 - (exp(-lam_foi[dn1[i]])) *
#                             exp(-sum(lam_inf[dn1[i]:(r[i] - 1)])))

#                 for (k in (dn1[i] + 1):(r[i] - 1)) {
#                     lik_temp <- lik_temp +
#                                 (1 - exp(-lam_foi[k])) *
#                                 exp(-sum(lam_sus[dn1[i]:(k - 1)])) *
#                                 exp(-sum(lam_inf[k:(r[i] - 1)])) *
#                                 exp(-sum(lam_foi[k:(r[i] - 1)]))
#                 }

#             } else {
#                 lik_temp <- (1 - exp(-lam_foi[dn1[i]])) *
#                             exp(-lam_inf[dn1[i]])
#             }
#             sumllik <- sumllik +
#                        log(exp(-sum(lam_sus[e[i]:(dn1[i] - 1)])) *
#                            exp(-sum(lam_foi[e[i]:(dn1[i] - 1)])) *
#                            lik_temp +
#                        exp(-sum(lam_foi[e[i]:(r[i] - 1)])) *
#                        exp(-sum(lam_sus[e[i]:(r[i] - 1)])))
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecNegCensPostNo = list(
#         BUGSdist = "dRecNegCensPostNo(n_samples,e,r,dn1,age2date,beta0_inf,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "dn1 = double(1)",
#             "age2date = double(1)",
#             "beta0_inf = double(0)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# assign("dRecNegCensPostNo", dRecNegCensPostNo, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dRecNegCensPostNo(
# #         x = 1,
# #         n_samples = n_fit_rec_neg_cens_postno,
# #         e = dat$d_fit_rec_neg_cens_postno$e_age,
# #         r = dat$d_fit_rec_neg_cens_postno$r_age,
# #         dn1 = dat$d_fit_rec_neg_cens_postno$rec_age,
# #         age2date = dat$d_fit_rec_neg_cens_postno$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end <- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected deer that were test neg at capture,
# ###   then test negative at recap,
# ###   that die
# ###
# ###   d_fit_rec_neg_mort
# ###
# ###   Overleaf Equation (18)
# ###
# #######################################################################

# dRecNegMort <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double(0)) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0
#             lam_sus <- 0

#             #############################################
#             ### preliminary hazards for the likelihood
#             #############################################
#             ### survival hazard for susceptible deer
#             for (j in e[i]:(r[i] - 1)) {
#                 lam_sus <- lam_sus + exp(beta0_sus +
#                     age_effect_surv[j] +
#                     period_effect_surv[(j + age2date[i])])
#             }
#             lam_sus_s <- sum(exp(beta0_sus +
#                 age_effect_surv[r[i]:(s[i] - 1)] +
#                 period_effect_surv[(r[i] + age2date[i]):
#                                    (s[i] - 1 + age2date[i])]))
#             ### force of infection infection hazard
#             for (j in e[i]:(s[i] - 1)) {
#                 lam_foi <- lam_foi + exp(foi_age_effect[age_lookup[j]])
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             sumllik <- sumllik - lam_sus + log(1 - exp(-lam_sus_s)) - lam_foi
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecNegMort = list(
#         BUGSdist = "dRecNegMort(n_samples,e,r,s,age2date,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "age2date =  double(1)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#               "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# # ###for a user-defined distribution
# assign("dRecNegMort", dRecNegMort, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dRecNegMort(x = 1,
# #         n_samples = n_fit_rec_neg_mort,
# #         e = dat$d_fit_rec_neg_mort$e_age,
# #         r = dat$d_fit_rec_neg_mort$r_age,
# #         s = dat$d_fit_rec_neg_mort$s_age,
# #         age2date = dat$d_fit_rec_neg_mort$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer that were test neg at capture,
# ###   then test positive at recap,
# ###   than these were right censored
# ###
# ###   d_fit_rec_pos_cens
# ###
# ###   Overleaf Equation (20)
# ###
# #######################################################################

# dRecPosCens <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0),
#                    e = double(1),
#                    r = double(1), 
#                    dn = double(1),
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {

#         sumllik <- 0 # intialize log-likelihood

#         for (i in 1:n_samples) {

#             lam_foi <- nimNumeric(r[i] - 1, init = FALSE)
#             lam_sus <- nimNumeric(r[i] - 1, init = FALSE)
#             lam_inf <- nimNumeric(r[i] - 1, init = FALSE)
#             lik_temp <- 0

#             #############################################
#             ### preliminary hazards for the likelihood
#             #############################################
#             ### survival hazard for susceptible deer
#             lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
#                 age_effect_surv[e[i]:(dn[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):(dn[i] - 1 + age2date[i])])

#             ### survival hazard while infected
#             lam_inf[e[i]:(r[i] - 1)] <- exp(beta0_inf +
#                 age_effect_surv[e[i]:(r[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (r[i] - 1 + age2date[i])])
#             ### force of infection infection hazard
#             lam_foi[e[i]:(dn[i] - 1)] <- exp(foi_age_effect[age_lookup[e[i]:(dn[i] - 1)]])

#             #######################################
#             ### calculating the joint likelihood
#             #######################################

#             lik_temp <- lik_temp + (1 - exp(-lam_foi[e[i]])) *
#                 exp(-sum(lam_inf[e[i]:(r[i] - 1)]))

#             for (k in (e[i] + 1):(dn[i] - 1)) {
#                 lik_temp <- lik_temp + 
#                     (1 - exp(-lam_foi[k])) *
#                     exp(-sum(lam_inf[k:(r[i] - 1)])) *
#                     exp(-sum(lam_foi[e[i]:(k - 1)])) *
#                     exp(-sum(lam_sus[e[i]:(k - 1)]))
#             }
#             sumllik <- sumllik + log(lik_temp)
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecPosCens = list(
#         BUGSdist = "dRecPosCens(n_samples,e,r,dn,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "dn = double(1)",
#             "age2date = double(1)",
#             "beta0_inf = double(0)",
#             "beta0_sus = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# # ### for a user-defined distribution
# assign("dRecPosCens", dRecPosCens, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dRecPosCens(
# #         x = 1,
# #         n_samples = n_fit_rec_pos_cens,
# #         e = dat$d_fit_rec_pos_cens$e_age,
# #         r = dat$d_fit_rec_pos_cens$r_age,
# #         dn = dat$d_fit_rec_pos_cens$rec_age,
# #         age2date =  dat$d_fit_rec_pos_cens$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   deer that were test neg at capture,
# ###   then test positive at recap,
# ###   than die
# ###
# ###   d_fit_rec_pos_mort
# ###
# ###   Overleaf Equation (22)
# ###
# #######################################################################

# dRecPosMort <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), 
#                    e = double(1), 
#                    r = double(1), 
#                    s = double(1), 
#                    dn = double(1), 
#                    age2date = double(1),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    foi_age_effect = double(1),
#                    age_lookup = double(1),
#                    log = double()) {

#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
#             lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
#             lik_temp <- 0

#             #############################################
#             ## preliminary hazards for the likelihood
#             #############################################
#             ### survival hazard for susceptible deer
#             lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
#                 age_effect_surv[e[i]:(dn[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                                     (dn[i] - 1 + age2date[i])])
#             ### survival hazard while infected
#             lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
#                 age_effect_surv[e[i]:(s[i] - 1)] +
#                 period_effect_surv[(e[i] + age2date[i]):
#                 (s[i] - 1 + age2date[i])])
#             ### force of infection infection hazard
#             lam_foi[e[i]:(dn[i] - 1)] <- 
#                 exp(foi_age_effect[age_lookup[e[i]:(dn[i] - 1)]])

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             lik_temp <- lik_temp +
#                         (1 - exp(-lam_foi[(e[i])])) *
#                 exp(-sum(lam_inf[(e[i]):(r[i] - 1)]))
#             if((dn[i] - e[i])>1){
#                 for (k in (e[i]+ 1):(dn[i] - 1)) {
#                     lik_temp <- lik_temp + 
#                         (1 - exp(-lam_foi[k])) *
#                         exp(-sum(lam_inf[k:(r[i] - 1)])) *
#                         exp(-sum(lam_foi[e[i]:(k - 1)])) *
#                         exp(-sum(lam_sus[e[i]:(k - 1)]))
#                 }
#             }
#             sumllik <- sumllik + 
#                        log(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) +
#                        log(lik_temp)
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecPosMort = list(
#         BUGSdist = "dRecPosMort(n_samples,e,r,s,dn,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "dn = double(1)",
#             "age2date = double(1)",
#             "beta0_sus = double(0)",
#             "beta0_inf = double(0)",
#             "age_effect_surv = double(1)",
#             "period_effect_surv = double(1)",
#             "foi_age_effect = double(1)",
#             "age_lookup = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ### Global Declaration so Nimble can access
# assign("dRecPosMort", dRecPosMort, envir = .GlobalEnv)


# # n_samples = n_fit_rec_pos_mort
# #         e = dat$d_fit_rec_pos_mort$e_age
# #         r = dat$d_fit_rec_pos_mort$r_age
# #         s = dat$d_fit_rec_pos_mort$s_age
# #         dn = dat$d_fit_rec_pos_mort$rec_age
# #         age2date = dat$d_fit_rec_pos_mort$age2date
# #         beta0_sus = beta0_survival_sus
# #         beta0_inf = beta0_survival_inf
# #         age_effect_surv = age_effect_true
# #         period_effect_surv = period_effect_true
# #         foi_age_effect = foi_age_effect
# #         age_lookup = age_lookup


# # starttime <- Sys.time()
# # test <-  dRecPosMort(
# #         x = 1,
# #         n_samples = n_fit_rec_pos_mort,
# #         e = dat$d_fit_rec_pos_mort$e_age,
# #         r = dat$d_fit_rec_pos_mort$r_age,
# #         s = dat$d_fit_rec_pos_mort$s_age,
# #         dn = dat$d_fit_rec_pos_mort$rec_age,
# #         age2date = dat$d_fit_rec_pos_mort$age2date,
# #         beta0_sus = beta0_survival_sus,
# #         beta0_inf = beta0_survival_inf,
# #         age_effect_surv = age_effect_true,
# #         period_effect_surv = period_effect_true,
# #         foi_age_effect = foi_age_effect,
# #         age_lookup = age_lookup,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test
