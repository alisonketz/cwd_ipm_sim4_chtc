#######################################
### Data for Model Fitting
#######################################

nimData <- list(
    Z_period = Z_period,
    Z_age = Z_age,
    y_fit_surv_neg = 1,
    fit_surv_neg_censored = d_fit_neg$rt_censor,
    fit_surv_neg_left = d_fit_neg$left_age,
    fit_surv_neg_right = d_fit_neg$right_age,
    fit_surv_neg_age2date = d_fit_neg$age2date,
    y_fit_surv_pos = 1,
    fit_surv_pos_censored = d_fit_pos$rt_censor,
    fit_surv_pos_left = d_fit_pos$left_age,
    fit_surv_pos_right = d_fit_pos$right_age,
    fit_surv_pos_age2date = d_fit_pos$age2date
    )

#######################################
### Constants for MCMC
#######################################

nimConsts <- list(
    n_year = n_year,
    n_ageclass = n_ageclass,
    nT_period = nT_period,
    nT_age = nT_age,
    n_yr_start_age = n_yr_start_age,
    nT_age_short = nT_age_short,
    nT_age_surv_aah = nT_age_surv_aah,
    n_age = n_age,
    age_lookup = age_lookup,
    yr_start_age = yr_start_age,
    n_yr_start_age = n_yr_start_age,
    yr_start_pop = yr_start_pop,
    nknots_age = nknots_age,
    nknots_period = nknots_period,
    n_fit_neg = n_fit_neg,
    n_fit_pos = n_fit_pos,
    # n_fit_sus_cens_posttest = n_fit_sus_cens_posttest,
    # n_fit_sus_cens_postno = n_fit_sus_cens_postno,
    # n_fit_sus_mort_posttest = n_fit_sus_mort_posttest,
    # n_fit_sus_mort_postno = n_fit_sus_mort_postno,
    # n_fit_icap_cens = n_fit_icap_cens,
    # n_fit_icap_mort = n_fit_icap_mort,
    # n_fit_rec_neg_cens_posttest = n_fit_rec_neg_cens_posttest,
    # n_fit_rec_neg_cens_postno = n_fit_rec_neg_cens_postno,
    # n_fit_rec_neg_mort = n_fit_rec_neg_mort,
    # n_fit_rec_pos_cens = n_fit_rec_pos_cens,
    # n_fit_rec_pos_mort = n_fit_rec_pos_mort,
    # n_fit_idead = n_fit_idead,
    nconst = 1 / sqrt(2 * pi)
)


#######################################
### Initial Values for MCMC
#######################################

initsFun <- function()list(
    beta0_sus_temp = rnorm(1, -5.75, 0.0001),
    sus_mix = 1,
    beta0_inf_temp = rnorm(1, -3, 0.0001),
    inf_mix = 1,
    ln_b_age_survival = rnorm(nknots_age) * 10^-4,
    tau_age_survival = runif(1, .1, 1),
    mix_survival = 1,
    ln_sk_period = runif(1, .1, 1),
    sda_period = runif(1, .1, 1),
    alpha_period = rnorm(nknots_period, 0, 1)#,
    # tau_age_foi = runif(1, 1.5, 1.7),
    # foi_age_effect = c(rnorm(1, -6.4, sd = .1),
    #               rnorm(1, -6.5, sd = .1),
    #               rnorm(1, -6.2, sd = .1),
    #               rnorm(1, -6.1, sd = .1),
    #               rnorm(1, -6, sd = .1),
    #               rnorm(1, -5, sd = .1))
    )
nimInits <- initsFun()

########################################################
### Build and compile model in R
########################################################

# start_Rmodel <- Sys.time()
Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun(),
                      calculate = FALSE,
                      check = FALSE
                      )
# end_Rmodel <- Sys.time() - start_Rmodel
Rmodel$initializeInfo()

Cnim <- compileNimble(Rmodel)

# for(i in 1:10){beepr::beep(1)}

#######################################
### Parameters to trace in MCMC
#######################################

parameters <- c(
            #   "tau_age_foi",
            #   "foi_age_effect",
            #   "foi_age_mu",
              "beta0_survival_inf",
              "beta0_survival_sus",
              "tau_age_survival",
              "age_effect_survival",
              "ln_b_age_survival",
              "mix_survival",
              "ln_sk_period",
              "sdk_period",
              "tauk_period",
              "stauk_period",
              "sda_period",
              "taua_period",
              "ratioinf_period",
              "period_effect_surv",
            #   "psi",
              "sn_inf",
              "sn_sus"
              # "sh_inf",
              # "sh_sus"#,
               )

confMCMC <- configureMCMC(Rmodel,
                         monitors = parameters,
                         thin = 1,
                         # enableWAIC = TRUE,
                         useConjugacy = FALSE)
nimMCMC <- buildMCMC(confMCMC)
CnimMCMC <- compileNimble(nimMCMC,
                          # resetFunctions = TRUE,
                         project = Rmodel)
# for(i in 1:10){beepr::beep(1)}

reps <- 10000
bin <- reps * .5
n_chains <- 3

# set.seed(1001)
starttime <- Sys.time()
mcmcout <- runMCMC(CnimMCMC,
                  niter = reps,
                  nburnin = bin,
                  nchains = n_chains,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )
runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")
runtime
# for (i in 1:10) {beepr::beep(1)}

# mcmcout$summary
# end_Rmodel
# endtime_rmodel_compile
# endtime_mcmc
# runtime

# sink("runtime_allsteps.txt")
# cat("Rmodel:\n")
# end_Rmodel
# cat("\nCompile Rmodel:\n")
# endtime_rmodel_compile
# cat("\nCompile MCMC:\n")
# endtime_mcmc
# cat("\nRun MCMC 100 iter: ",runtime)
# sink()
#############################################################
###
### Running in parallel, restartable
###
#############################################################

# reps  <- 10
# bin <- reps * .5
# n_thin <- 1
# n_chains <- 3
# starttime <- Sys.time()
# cl <- makeCluster(n_chains, timeout = 5184000)

# clusterExport(cl, c("modelcode",
#                     "initsFun",
#                     "nimData",
#                     "nimConsts",
#                     "parameters",
#                     "reps",
#                     "bin",
#                     "n_thin",
#                     "set_period_effects_constant",
#                     "dInfHarvest",
#                     "dSusHarvest",
#                     "dSusCensTest",
#                     "dSusCensNo",
#                     "dSusMortTest",
#                     "dSusMortNoTest",
#                     "dIcapCens",
#                     "dIcapMort",
#                     "dRecNegCensTest",
#                     "dRecNegMort",
#                     "dRecPosMort",
#                     "dRecPosCens",
#                     "dNegCapPosMort",
#                     "dAAH",
#                     "calc_surv_aah",
#                     "calc_surv_harvest",
#                     "calc_infect_prob"
#                     ))
# for (j in seq_along(cl)) {
#   set.seed(j + 1000)
#   init <- initsFun()
#   clusterExport(cl[j], "init")
# }
# for (i in 1:10) {beepr::beep(1)}

# # starttime <- Sys.time()
# # mcmcout1 <-  mcmc.list(clusterEvalQ(cl, {
# #   library(nimble)
# #   library(coda)

#   compile(dInfHarvest)
#   compile(dSusHarvest)
#   compile(dSusCensTest)
#   compile(dSusCensNo)
#   compile(dSusMortTest)
#   compile(dSusMortNoTest)
#   compile(dIcapCens)
#   compile(dIcapMort)
#   compile(dRecNegCensTest)
#   compile(dRecNegMort)
#   compile(dRecPosMort)
#   compile(dRecPosCens)
#   compile(dNegCapPosMort)
#   compile(dAAH)
#   compile(set_period_effects_constant)
#   compile(calc_surv_aah)
#   compile(calc_surv_harvest)
#   compile(calc_infect_prob)

# #   ##############################################################
# #   ###
# #   ### Execute MCMC
# #   ###
# #   ##############################################################

# #   Rmodel <- nimbleModel(code = modelcode,
# #                         name = "modelcode",
# #                         constants = nimConsts,
# #                         data = nimData,
# #                         inits = initsFun)
# #   Cnim <- compileNimble(Rmodel)
# #   confMCMC <- configureMCMC(Rmodel,
# #                             monitors = parameters,
# #                             thin = n_thin,
# #                             useConjugacy = FALSE)
# #   nimMCMC <- buildMCMC(confMCMC)
# #   CnimMCMC <- compileNimble(nimMCMC,
# #                             project = Rmodel)

# #   CnimMCMC$run(reps, reset = FALSE)

# #   return(as.mcmc(as.matrix(CnimMCMC$mvSamples)))
# # }))
# # runtime1 <- difftime(Sys.time(), starttime, units = "min")

# # for(chn in 1:nc) { # nc must be > 1
# #   ind.keep <- c()
# #   for(p in 1:length(parameters)) ind.keep <-
# #       c(ind.keep, which(str_detect(dimnames(out1[[chn]])[[2]], parameters[p]))) %>% unique()
# #   out1[[chn]] <- out1[[chn]][,ind.keep]
# # }

# # ## Check convergence ##
# # out2 <- out1
# # ni.saved <- nrow(out2[[1]])
# # for(chn in 1:nc) { # nc must be > 1
  
# #   if(nb < 1) {
# #     nb.real <- (round(ni.saved * nb)+1)
# #   } else {
# #     nb.real <- (round(nb/nt)+1)
# #   }
# #   out2[[chn]] <- out2[[chn]][nb.real:ni.saved,]
# # }
# # out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))
# # stopCluster(cl)

# # save(mcmcout1, file = "mcmcout1.Rdata")
# # save(runtime1, file = "runtime1.Rdata")
# # save(endtime_rmodel_compile, file = "endtime_rmodel_compile.Rdata")
# # save(endtime_mcmc, file = "endtime_mcmc.Rdata")

# #not calculating waic, because too many params would need to be traced
# # posteriorSamplesMatrix <- rbind(mcmcout[[1]], mcmcout[[2]], mcmcout[[3]])
# # CnimMCMC$run(5)   ## non-zero number of iterations
# # nimble:::matrix2mv(posteriorSamplesMatrix, CnimMCMC$mvSamples)
# # # CnimMCMC$enableWAIC <- TRUE
# # waic_spline <- calculateWAIC(posteriorSamplesMatrix,Rmodel)


# # waic_spline_covs <- mcmcout$WAIC
# # save(waic_spline, file = "waic_spline.Rdata")



# ###
# ### save model run
# ###

# # save(runtime,file="results/runtime.Rdata")
# # save(mcmcout,file="results/mcmcout.Rdata")