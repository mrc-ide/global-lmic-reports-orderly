#https://doi.org/10.1016/S1473-3099(22)00001-9
#https://github.com/will-s-hart/variant_generation_times
#download matlab file from repo
download.file("https://github.com/will-s-hart/variant_generation_times/raw/main/Results/mcmc_posterior_mech.mat", "posteriors/temp", mode = "wb")
posterior <-  R.matlab::readMat("posteriors/temp")
unlink("posteriors/temp")
#extract the values we need
delta_gen_time_mean <- posterior$mean.post2
delta_gen_time_sd <- posterior$sd.post2
saveRDS(
  tibble(mean = delta_gen_time_mean, sd = delta_gen_time_sd),
  "posteriors/hart_et_al_posterior.Rds"
)
