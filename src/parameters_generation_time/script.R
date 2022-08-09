library(tidyverse)
library(drjacoby)
odin_infection_model <- odin::odin({
  initial(E1) <- 1
  initial(E2) <- 0
  initial(ICase1) <- 0
  initial(ICase2) <- 0
  initial(IMild) <- 0

  deriv(E1) <- -E1*gamma_E
  deriv(E2) <- (E1 - E2)*gamma_E
  deriv(ICase1) <- E2 * gamma_E - ICase1 * gamma_ICase
  deriv(ICase2) <- (ICase1 - ICase2) * gamma_ICase
  deriv(IMild) <- E2 * gamma_E - IMild * gamma_IMild

  gamma_E <- user()
  gamma_ICase <- user()
  gamma_IMild <- user()
  prop_Mild <- user()

  output(infectious) <- (ICase1 + ICase2) * (1 - prop_Mild) + prop_Mild * IMild
})
#parameters basics + multiplers for each variant + proportion mild for each study
Rcpp::sourceCpp("funcs.cpp")
##data, have to generate our data from studies etc
n_data <- 100*50
quants <- seq(0, 1, length.out = n_data + 2)[-c(1, n_data + 2)]
#wild type
#https://doi.org/10.1016/S1473-3099(20)30287-5
#gives distribution for serial interval, have to assume its the same as the gen time
#we'll use samples for this this
bi_data <- qgamma(quants, 2.29, 0.36)
#https://doi.org/10.1101/2020.09.04.20188516
#gives distribution for the generation time
shape <- (1.8/5.5)^(-1.086)
scale <- 5.5/gamma(1 + 1/shape)
ferreti_data <- qweibull(quants, shape, scale)
##Other, either use to make synthetic data or as quasi-prior distributions
#Wild
#10.2807/1560-7917.ES.2020.25.17.2000257
#gives a distribution for the parameters of the GI given that its Gamma dist
#use these to generate 100 samples from  thinned sampls
ganyani_posterior <- readRDS("posteriors/ganyani_et_al_posterior.Rds")
ganyani_posterior <- ganyani_posterior[round(seq(1, nrow(ganyani_posterior), length.out = n_data/100)),]
ganyani_data <- unlist(map(transpose(ganyani_posterior), function(pars){
  rate <- pars$mean/pars$sd
  qgamma(seq(0, 1, length.out = 102)[-c(1, 102)], rate = rate, shape = pars$mean * rate)
}))
#Delta
#https://doi.org/10.1016/S1473-3099(22)00001-9
#gives distributions for mean and sd of generation times
#use the samples to produce a density estimate of the underlying distributions and
#fit to that likelihood
hart_posterior <- readRDS("posteriors/hart_et_al_posterior.Rds")
hart_posterior <- hart_posterior[round(seq(1, nrow(hart_posterior), length.out = n_data/100)),]
hart_data <- unlist(map(transpose(hart_posterior), function(pars){
  #assume gamma
  rate <- pars$mean/pars$sd
  qgamma(seq(0, 1, length.out = 102)[-c(1, 102)], rate = rate, shape = pars$mean * rate)
}))
#Omicron
#https://doi.org/10.1016/j.lanepe.2022.100446
#gamma distributed, gives distributions for parameters
#having trouble getting their C code to compile, use central for now
manica_data <- qgamma(quants, shape = 2.39, scale = 2.95)
#https://doi.org/10.3390/v14030533
#normal distributed? serial interval + no code, leave it out

##Setup objects
params <- drjacoby::define_params(
    name = "dur_E", min = 0.5, max = 15, block = 1,
    name = "dur_ICase", min = 0.5, max = 15, block = 1,
    name = "dur_IMild", min = 0.5, max = 15, block = 1,
    name = "dur_E_delta", min = 0.5, max = 15, block = 2,
    name = "dur_ICase_delta", min = 0.5, max = 15, block = 2,
    name = "dur_IMild_delta", min = 0.5, max = 15, block = 2,
    name = "dur_E_omicron", min = 0.5, max = 15, block = 3,
    name = "dur_ICase_omicron", min = 0.5, max = 15, block = 3,
    name = "dur_IMild_omicron", min = 0.5, max = 15, block = 3
)
data <- list(wild = sort(c(bi_data, ferreti_data, ganyani_data)),
             delta = sort(c(hart_data)),
             omicron = sort(c(manica_data)))
#precalculate these for speed
t_max <- 100
N <- 100
denom_values <-  seq(0, t_max, length.out = N)
trapezoid_multipler <- log(t_max/N/2)
misc <- list(
  model_instance = odin_infection_model$new(user = list(
    prop_Mild = 2, gamma_E = 2/10, gamma_ICase = 2/10, gamma_IMild = 1/5
  )),
  denom_values = denom_values,
  trapezoid_multipler = trapezoid_multipler
)

#run the mcmc
mcmc_out <- drjacoby::run_mcmc(
  data = data, df_params = params, misc = misc, loglike = "loglike",
  logprior = "logprior", burnin = 1000, samples = 1000, rungs = 50, chains = 5
)
drjacoby::plot_par(mcmc_out)
drjacoby::plot_mc_acceptance(mcmc_out)
