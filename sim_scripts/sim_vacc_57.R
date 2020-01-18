prop_virus_sens <- 0.833
ve_24_months    <- 0.6
prop_coverage   <- 0.9
risk_comp_rr    <- 0.7
model_name      <- "vacc_sens_0.833_ve_0.6_cov_0.9_rc_0.7"

saved_dat <- paste0("/gscratch/csde/kpeebles/dat_prop_sens_", prop_virus_sens, "_cov_", prop_coverage, ".RDATA")

## Attach packages
library(evonet)
library(EpiModelHPC)

## Load default parameters
primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

## Combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)

## Restart parameters
evoparams$restart_val  <- "restart_vacc_cea"
evoparams$restart_time <- 365 * 29

## Override selected default parameters
evoparams$initial_pop            <- 30000
evoparams$initial_infected       <- 60
evoparams$model_sex              <- "hetero"
evoparams$n_steps                <- 365 * 44
evoparams$popsumm_frequency      <- 10
evoparams$poisson_birth_lambda   <- 0.0137*(evoparams$initial_pop/100)
evoparams$vl_peak_agent_flag     <- T
evoparams$fast_edgelist          <- T
evoparams$mean_sex_acts_day      <- 0.36
evoparams$min_age                <- 16
evoparams$max_age                <- 55
evoparams$sti_prob               <- 0.692
evoparams$asmr_data_male         <- "south_africa_male_1990"
evoparams$asmr_data_female       <- "south_africa_female_1990"
evoparams$initial_agedata_male   <- "south_africa_male_16_to_100_1990"
evoparams$initial_agedata_female <- "south_africa_female_16_to_100_1990"
evoparams$trans_RR_STI           <- 2.14
evoparams$AverageLogSP0          <- 4.83
evoparams$trans_lambda           <- 0.000720
evoparams$prop_AI                <- 0.077
evoparams$mean_prop_acts_AI      <- 0.235
evoparams$sd_prop_acts_AI        <- 0.120

evoparams$nsims                  <- 64
evoparams$ncores                 <- 16
evoparams$output_path            <- '/gscratch/csde/kpeebles/vacc_sims'

## Circumcision parameters
evoparams$circum_prob        <- 0.3
evoparams$circum_prob_chg    <- c(0.45,   0.5,    0.55,   0.6,    0.65,   0.7,    0.75,   0.85)
evoparams$circum_prob_yr_chg <- c(12*365, 13*365, 15*365, 16*365, 18*365, 20*365, 22*365, 24*365)

## Condom use parameters
evoparams$condom_prob        <- 0 # This value will be used only for coital acts at first timestep. Value will be reassigned by Hill function in subsequent timesteps.
evoparams$condom_prob_change <- T

## Risk compensation in form of decreased condom use among those who are vaccinated
evoparams$risk_comp_cond    <- T
evoparams$risk_comp_cond_rr <- risk_comp_rr

## Treatment parameters
evoparams$trtmnt_sex_age            <- T
evoparams$cd4_treatment_threshold   <- 4 # Treatment eligibility is initially <200
evoparams$start_treatment_campaign  <- c(0, 21*365, 24*365, 26*365) # Treatment eligibility changes at years 21, 24, and 26
evoparams$cd4_trt_guidelines_chgs   <- list(4, 3:4, 2:4, 1:4) # Treatment eligibility changes over time from <200, <350, <500, and all +
evoparams$tx_in_acute_phase         <- T
evoparams$cov_prob                  <- c(0, 0.010, 0.021, 0.030, 0.049, 0.100, 0.191, 0.283, 0.402, 0.560) # 0.56 from UNAIDS estimate, 2016
evoparams$cov_prob_yrs              <- c(0, 11:18, 27) # Years at which coverage changes 
evoparams$cov_prob_scal <- matrix(c(0.569, 1.240, 1.240, 0.421, 0.919, 0.919), ncol = 2, 
                                  dimnames = list(c("15-24", "25-34", "35+"), c("f", "m")))
evoparams$cov_prob_ageg <- list(c(15, 25), c(25, 35), c(35, evoparams$max_age + 1))

## Vaccination parameters
evoparams$vacc_wane               <- T
evoparams$start_vacc_campaign     <- (29*365):(44*365)
evoparams$perc_virus_vaccine_sens <- prop_virus_sens # prop_virus_sens * 100% of virus is sensitive and VE at 24 months is ve_24_months
evoparams$ve_24_months            <- ve_24_months # Value used in waning vaccine-induced immunity function to determine daily RR against sensitive virus
evoparams$perc_vaccinated         <- 1/(3*365) # daily probability of vaccination is equal to prob_care/time to achieve target coverage
evoparams$vacc_eff_duration       <- 365*5 # After 5 years, vaccinated individuals will receive a booster vaccine and immunity will return to initial value
evoparams$prob_care               <- prop_coverage # eligible_care == 1 is vaccination criteria, so a maximum of prop_coverage * 100% susceptible agents will be vaccinated

## Network formation parameters
evoparams$nw_form_terms <- "~edges + concurrent('sex') + absdiffby('age', 'sex', 3, c('f', 'm')) + offset(nodematch('sex', diff=FALSE))"
evoparams$nw_coef_form  <- -Inf
nEdges                  <- 0.76*(evoparams$initial_pop/2)
evoparams$target_stats  <- c(nEdges, 0.079*(evoparams$initial_pop/2), 0.176*(evoparams$initial_pop/2), nEdges)
evoparams$relation_dur  <- 1409

## Add parameters that are functions of other input parameters
evoparams          <- input_parameters_derived(evoparams)
evoparams$age_dist <- seq(50, 10, -10/9)/1110

## Convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)

## Initialize network
nw <- setup_initialize_network(evoparams)

## Run QAQC on input parameters
input_parameters_qaqc(evoparams)

## Estimate initial network model
netest_arg_list <- list(
  nw            =  nw,
  formation     =  as.formula(evoparams$nw_form_terms),
  target.stats  =  evoparams$target_stats,
  coef.form     =  evoparams$nw_coef_form,
  constraints   =  as.formula(evoparams$nw_constraints),
  verbose       =  FALSE,
  coef.diss     =  dissolution_coefs(dissolution =  as.formula(evoparams$dissolution),
                                     duration    =  evoparams$relation_dur,
                                     d.rate      =  3e-05),
  set.control.ergm = control.ergm(MCMLE.maxit = 100))

estimated_nw <- do.call(EpiModel::netest, netest_arg_list)

## Create initial vector of infection status as an EpiModel object
infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                    status.rand = FALSE)

## Create list of modules to run for input into epimodel_control_fxn() below
evo_module_list <- list(
  "initialize.FUN"     = initialize_module,
  "restart.FUN"        = restart_module,
  "aging.FUN"          = aging,
  "treatment.FUN"      = social_treatment_sex_age,
  "vacc.FUN"           = social_treatment_vaccination_waning,
  "update_vl.FUN"      = viral_update_gamma,
  "update_cd4.FUN"     = viral_update_cd4_daily,
  "coital_acts.FUN"    = social_coital_acts_module,
  "trans.FUN"          = transmission_main_module,
  "trans_book.FUN"     = transmission_bookkeeping_module,
  "trans_cd4.FUN"      = transmission_cd4_module,
  "deaths.FUN"         = vital_deaths_module,
  "births.FUN"         = vital_births_module,
  "summary.FUN"        = summary_module,
  "resim_nets.FUN"     = EpiModel::resim_nets,
  "verbose.FUN"        = NULL)

## Call EpiModel's control fxn (load evonet modules into EpiModel)
evocontrol <- setup_epimodel_control_object(evonet_params = evoparams,
                                            module_list   = evo_module_list)

## Run simulation
evomodel <- EpiModel::netsim(x = estimated_nw,
                             param = evoparams,
                             init = infected_list,
                             control = evocontrol)

plots_popsumm(evomodel, outpath = evoparams$output_path,
              name = model_name, nw_stats = TRUE, max_points_rep = 100,
              evoparams$popsumm_frequency)

evomodel <- pop_vars(model = evomodel)

## Save simulation
evomodel$epi              <- NULL
evomodel$stats            <- NULL
evomodel$control          <- NULL
evomodel$attr             <- NULL
evomodel$nw               <- NULL
evomodel$coital_acts_list <- NULL
evomodel$vl_list          <- NULL
evomodel$InfMat           <- NULL
evomodel$age_list         <- NULL
evomodel$el               <- NULL
evomodel$partner_list     <- NULL

assign(model_name, evomodel)

file_name <- paste(model_name,".RData",sep="")

save(list = model_name,
     file = file.path(evoparams$output_path, file_name))
