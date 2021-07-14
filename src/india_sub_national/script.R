RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

system(paste0("echo India Subnational for  ",state))

## -----------------------------------------------------------------------------
## 0. Checks and Function Definitions
## -----------------------------------------------------------------------------

# check on the state name
if (!(state %in% c(
  "Andaman and Nicobar Islands","Andhra Pradesh","Arunachal Pradesh","Assam",
  "Bihar","Chandigarh","Chhattisgarh","Dadra and Nagar Haveli and Daman and Diu",
  "Delhi","Goa","Gujarat","Haryana","Himachal Pradesh","Jammu and Kashmir",
  "Jharkhand","Karnataka","Kerala","Ladakh","Lakshadweep","Madhya Pradesh",
  "Maharashtra","Manipur","Meghalaya","Mizoram","Nagaland","Odisha","Puducherry",
  "Punjab","Rajasthan","Sikkim","Tamil Nadu","Telangana","Tripura","Uttar Pradesh",
  "Uttarakhand","West Bengal"))) {
  stop("State is not correct")
}

## -----------------------------------------------------------------------------
## 1. GET INPUT DATA
## -----------------------------------------------------------------------------

## a. Get from local files
## -----------------------------------------------------------------------------

# get pop data from file
demog <- readRDS("demog.rds")
pop <- demog$n[demog$state == state]

# get icu beds from file
icu_beds <- readRDS("icu_beds.rds")
icu_beds <- icu_beds$icu_beds[icu_beds$state == state]

# get hosp beds from file
hosp_beds <- readRDS("hosp_beds.rds")
hosp_beds <- hosp_beds$hosp_beds[hosp_beds$state == state]

sero_df <- readRDS("sero.rds")
sero_df <- sero_df[sero_df$state == state,]

## b. Get from remote changing sources
## -----------------------------------------------------------------------------

# get death data
subnat_df <- read.csv("https://api.covid19india.org/csv/latest/states.csv") %>% 
  filter(!(State %in% c("India", "State Unassigned"))) %>% 
  mutate(Date = as.Date(Date)) %>% 
  group_by(State) %>% 
  complete(Date = seq.Date(min(as.Date(Date)), max(as.Date(Date)), 1)) %>% 
  mutate(State = replace_na(State, unique(!is.na(State)))) %>% 
  mutate(cases = replace_na(Confirmed, 0),
         deaths = replace_na(Deceased, 0)) %>% 
  select(-c("Tested", "Recovered", "Other", "Confirmed", "Deceased")) %>% 
  rename(date = Date, state = State)
df <- subnat_df[subnat_df$state == state, ] %>% 
  ungroup %>% 
  select(date, deaths, cases) %>% 
  arrange(date)
df$deaths <- c(df$deaths[1], diff(df$deaths))
df$deaths[df$deaths < 0] <- 0
df$cases <- c(df$cases[1], diff(df$cases))
df$cases[df$cases < 0] <- 0

# filter to the current date
df <- df[as.Date(df$date) <= as.Date(date),]

# get vaccination data
subnat_vacc <- read.csv("https://raw.githubusercontent.com/sociepy/covid19-vaccination-subnational/main/data/countries/India.csv")
subnat_vacc <- subnat_vacc %>% 
  mutate(region = replace(region, which(region_iso %in% c("IN-DN", "IN-DD")), "Dadra and Nagar Haveli and Daman and Diu"),
         region_iso = replace(region_iso, which(region_iso %in% c("IN-DN", "IN-DD")), "IN-DD")) %>% 
  group_by(region, date, region_iso) %>% 
  summarise(total_vaccinations = sum(total_vaccinations),
            people_vaccinated = sum(people_vaccinated),
            people_fully_vaccinated = sum(people_fully_vaccinated)) %>% 
  rename(state = region)
subnat_vacc <- subnat_vacc[subnat_vacc$state == state, ]
vacc_inputs <- get_vaccine_inputs(max(df$date), subnat_vacc)

## -----------------------------------------------------------------------------
## 2. Fit Model
## -----------------------------------------------------------------------------

# fit model
res <- fit_spline_rt(data = df, 
                     country = as.character("India"), 
                     pop = pop, 
                     reporting_fraction = as.numeric(rf), 
                     vacc_inputs = vacc_inputs,
                     n_mcmc = as.numeric(n_mcmc),
                     replicates = as.numeric(replicates),
                     hosp_beds = as.numeric(hosp_beds),
                     icu_beds = as.numeric(icu_beds),
                     sero_df = sero_df) 


# add state for ease and remove the output for memory
res$parameters$state <- state
output <- res$output
res$output <- NULL

# save output without output for memory
saveRDS(res, "res.rds")
res$output <- output

# make a quick plot so we can check fits easily
rtp <- rt_plot_immunity(res)
dp <- plot(res, particle_fit = TRUE) + theme(legend.position = "none")
sero <- sero_plot(res, sero_df)
ar <- ar_plot(res)

rf_over <- paste0(round(quantile(res$replicate_parameters$rf)[c(2,4)], digits = 2)*100, "%", collapse = " - ")
ggsave("fitting.pdf",width=12, height=12, 
       cowplot::plot_grid(rtp$plot + ggtitle(paste0(state, ". Death Reporting at ", rf_over)), 
                          dp, sero, ar, ncol = 1))


