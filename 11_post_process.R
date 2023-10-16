###############################################################################################
###
### Post processing
###
###############################################################################################

fit_sum <- mcmcout$summary$all.chains
out <- mcmcout$samples

modelid <- "A"

#####################
###
### Objects to save
###
#####################

true_values <- list(
     dat = dat,
     psi_true = psi_true,
     sn_sus_true = sn_sus_true,
     sn_inf_true = sn_inf_true
)

gd <- gelman.diag(out[, c(grep("beta", rownames(fit_sum)),
                          grep("sd", rownames(fit_sum)),
                          grep("foi_age_effect", rownames(fit_sum)),
                          grep("tau", rownames(fit_sum)))
                          ], multivariate = FALSE)

# waic <- mcmcout$WAIC

ess <- effectiveSize(out[, c(grep("beta", rownames(fit_sum)),
                          grep("sd", rownames(fit_sum)),
                          grep("foi_age_effect", rownames(fit_sum)),
                          grep("tau", rownames(fit_sum)))
                          ])


###############################################
###
### Save results
###
##############################################

save(fit_sum, file = paste0("fit_sum_",modelid,"_", processnum + 1,".Rdata"))
save(gd, file = paste0("gd_",modelid,"_", processnum + 1, ".Rdata"))
save(runtime, file = paste0("runtime_",modelid,"_", processnum + 1, ".Rdata"))
save(ess, file = paste0("ess_",modelid,"_", processnum + 1, ".Rdata"))
save(reps, file = paste0("reps_",modelid,"_", processnum + 1, ".Rdata"))
save(seedy, file = paste0("seed_",modelid,"_", processnum + 1, ".Rdata"))
save(true_values, file = paste0("true_values_",modelid,"_", processnum + 1, ".Rdata"))

# save(fit_sum, file = paste0("results/",modelid,"/fit_sum_", processnum + 1, ".Rdata"))
# save(gd, file = paste0("results/",modelid,"/gd_", processnum + 1, ".Rdata"))
# save(runtime, file = paste0("results/",modelid,"/runtime_", processnum + 1, ".Rdata"))
# save(ess, file = paste0("results/",modelid,"/ess_", processnum + 1, ".Rdata"))
# save(reps, file = paste0("results/",modelid,"/reps_", processnum + 1, ".Rdata"))
# save(seedy, file = paste0("results/",modelid,"/seed_", processnum + 1, ".Rdata"))

# save(true_values, file = paste0("true.values_", processnum + 1, ".Rdata"))
# save(mcmcout,file=paste0("mcmcout_",processnum+1,".Rdata"))
# save(waic, file = paste0("waic_", processnum + 1, ".Rdata"))
# save(agetime.out, file = paste0("agetime.out_", processnum + 1, ".Rdata"))

#############################
### Saving Model Description
#############################

# sink(paste0("figures/", modelid, "/model_description_", modelid, ".txt"))
sink(paste0("model_description_", modelid,'_', processnum + 1, ".txt"))
cat("Model description specifics:\n")
cat("niter:  ", reps, "\n")
cat("burnin:  ", bin, "\n")
cat("n_chains:  ", n_chains, "\n")
cat('seed: "', seedy, "\n")
cat("Model Variation:\n")
cat("no FOI period effects\n")
cat("includes FOI age effects RW2 with implicit interept \n\n")
cat("includes survival period effects kernel convolution\n")
cat("includes survival age effects cgam convex \n\n")
cat("runtime:  ", runtime, "\n")
cat("Gelman diag:")
print(gelman.diag(out[, grep("beta", rownames(fit_sum))], multivariate = FALSE))
print(gelman.diag(out[, grep("tau", rownames(fit_sum))], multivariate = FALSE))
print(gelman.diag(out[, grep("sd", rownames(fit_sum))], multivariate = FALSE))
# print(gelman.diag(out[, grep("foi_age_effect", rownames(fit_sum))], multivariate = FALSE))
cat("\nSummary Stats:  \n")
print(fit_sum)
sink()
