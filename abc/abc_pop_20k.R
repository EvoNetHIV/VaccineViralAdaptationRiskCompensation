## Script to run approximate Bayesian computation calibration for the South African heterosexual model

## Attach packages
library(evonet)
library(EasyABC)
library(EpiModelHPC)

## Function to be passed to ABC_sequential function
get_params <- function(x) {
  set.seed(x[1])
  require(evonet)
  
  ## Load default parameters
  primary_parameters  <- input_parameters_primary()
  cd4_data            <- input_parameters_cd4_data()
  
  ## Combine individual parameters into single list
  evoparams <- c(primary_parameters, cd4_data)
  
  ## Override selected default parameters
  evoparams$initial_pop            <- 20000
  evoparams$initial_infected       <- 40
  evoparams$model_sex              <- "hetero"
  evoparams$n_steps                <- 365*25
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
  evoparams$output_path            <- '/gscratch/csde/kpeebles'
  evoparams$trans_RR_STI           <- 2.14
  evoparams$AverageLogSP0          <- 4.83
  evoparams$trans_lambda           <- x[2]
  evoparams$prop_AI                <- x[3]
  evoparams$mean_prop_acts_AI      <- x[4]
  evoparams$sd_prop_acts_AI        <- x[5]
  
  ## Circumcision parameters
  evoparams$circum_prob            <- 0.3
  evoparams$circum_prob_chg        <- c(0.45,   0.5,    0.55,   0.6,    0.65,   0.7,    0.75,   0.85)
  evoparams$circum_prob_yr_chg     <- c(12*365, 13*365, 15*365, 16*365, 18*365, 20*365, 22*365, 24*365)
  
  ## Condom use parameters
  evoparams$condom_prob            <- 0 # This value will be used only for coital acts at first timestep. Value will be reassigned by Hill function in subsequent timesteps.
  evoparams$condom_prob_change     <- T
  
  ## Treatment parameters
  evoparams$trtmnt_sex_age            <- T
  evoparams$cd4_treatment_threshold   <- 4 # Treatment eligibility is initially <200
  evoparams$start_treatment_campaign  <- c(0, 21*365, 24*365, 26*365) # Treatment eligibility changes at years 21, 24, and 26
  evoparams$cd4_trt_guidelines_chgs   <- list(4, 3:4, 2:4, 1:4) # Treatment eligibility changes over time from <200, <350, <500, and all +
  evoparams$tx_in_acute_phase         <- T
  evoparams$cov_prob                  <- c(0, 0.01, 0.021, 0.030, 0.049, 0.100, 0.191, 0.283, 0.402, 0.78)
  evoparams$cov_prob_yrs              <- c(0, 11:18, 23) # Years at which coverage changes 
  evoparams$cov_prob_scal <- matrix(c(0.718571, 1.000558, 1.239242, 0.584926, 0.814467, 1.008758), ncol = 2, 
                                    dimnames = list(c("15-24", "25-34", "35+"), c("f", "m")))
  evoparams$cov_prob_ageg <- list(c(15, 25), c(25, 35), c(35, evoparams$max_age + 1))
  
  ## Network formation parameters
  evoparams$nw_form_terms <- "~edges + concurrent('sex') + absdiffby('age', 'sex', 3, c('f', 'm')) + offset(nodematch('sex', diff = F))"
  evoparams$nw_coef_form  <- -Inf
  nEdges                  <- 0.76*(evoparams$initial_pop/2)
  fem_conc                <- x[6]
  male_conc               <- x[7]
  evoparams$target_stats  <- c(nEdges, fem_conc*(evoparams$initial_pop/2), male_conc*(evoparams$initial_pop/2), nEdges)
  evoparams$relation_dur  <- x[8]

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
                                       d.rate      =  3e-05))
  
  estimated_nw <- do.call(EpiModel::netest, netest_arg_list)
  
  ## Create initial vector of infection status as an EpiModel object
  infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                      status.rand = FALSE)
  
  ## Create list of modules to run for input into epimodel_control_fxn() below
  evo_module_list <- list(
    "initialize.FUN"     = initialize_module,
    "aging.FUN"          = vital_aging_module,
    "treatment.FUN"      = social_treatment_sex_age,
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
  
  ## Calculate prevalence summary statistic from simulated data at years for which empirical data is available
  sa.year   <- c(seq(1,13,1), 16, 19, 23) # Year 1 corresponds to end of 1990, 16 to end of 2005
  model.year <- sa.year * 365/evomodel$param[[1]]$popsumm_freq
  
  out1 <- evomodel$popsumm[[1]]$prev_15to49[model.year] # Ages 15-49, 1990-2001, 2002, 2005, 2008, 2012
  
  out2 <- c(evomodel$popsumm[[1]]$prev_m_15to24[13 * 365/evomodel$param[[1]]$popsumm_freq], # 2002, males 15-24
            evomodel$popsumm[[1]]$prev_f_15to24[13 * 365/evomodel$param[[1]]$popsumm_freq], # 2002, females 15-24
            evomodel$popsumm[[1]]$prev_m_15to49[13 * 365/evomodel$param[[1]]$popsumm_freq], # 2002, males 15-49
            evomodel$popsumm[[1]]$prev_f_15to49[13 * 365/evomodel$param[[1]]$popsumm_freq]) # 2002, females 15-49
  
  out3 <- c(evomodel$popsumm[[1]]$prev_m_15to24[16 * 365/evomodel$param[[1]]$popsumm_freq], # 2005, males 15-24
            evomodel$popsumm[[1]]$prev_f_15to24[16 * 365/evomodel$param[[1]]$popsumm_freq], # 2005, females 15-24
            evomodel$popsumm[[1]]$prev_m_15to49[16 * 365/evomodel$param[[1]]$popsumm_freq], # 2005, males 15-49
            evomodel$popsumm[[1]]$prev_f_15to49[16 * 365/evomodel$param[[1]]$popsumm_freq]) # 2005, females 15-49
  
  out4 <- evomodel$popsumm[[1]]$prev_15to24[19 * 365/evomodel$param[[1]]$popsumm_freq] # 2008, 15-24
  
  out5 <- c(evomodel$popsumm[[1]]$prev_m_15to24[23 * 365/evomodel$param[[1]]$popsumm_freq], # 2012, males 15-24
            evomodel$popsumm[[1]]$prev_f_15to24[23 * 365/evomodel$param[[1]]$popsumm_freq], # 2012, females 15-24
            evomodel$popsumm[[1]]$prev_m_15to49[23 * 365/evomodel$param[[1]]$popsumm_freq], # 2012, males 15-49
            evomodel$popsumm[[1]]$prev_f_15to49[23 * 365/evomodel$param[[1]]$popsumm_freq]) # 2012, females 15-49
  
  out <- c(out1, out2, out3, out4, out5)
  
  return(out)
}

## Specify priors for per-act infectivity and relationship duration
priors  <- list(c("normal", 0.0005005, 0.00015), # trans_lambda
                c("unif", 0.01, 0.15),           # prop_AI 
                c("unif", 0.01, 0.60),           # mean_prop_acts_AI
                c("unif", 0.05, 0.30),           # sd_prop_acts_AI
                c("unif", 0.01, 0.15),           # fem_conc
                c("unif", 0.07, 0.20),           # male_conc
                c("unif", 254, 1971))            # rel_dur

## Specify prevalence targets for ABC fitting procedure
sa.prev <- c(.002, .005, .010, .019, .031, .048, .067, .088, .108, .126, .141, .153, .156, .162, .169, .188, # Ages 15-49, 1990-2001, 2002, 2005, 2008, 2012
             0.061, 0.120, 0.128, 0.177, # 2002: males, 15-24; females, 15-24; males, 15-49; females, 15-49
             0.044, 0.169, 0.117, 0.202, # 2005: males, 15-24; females, 15-24; males, 15-49; females, 15-49
             0.087,                      # 2008: 15-24
             0.029, 0.114, 0.145, 0.232  # 2012: males, 15-24; females, 15-24; males, 15-49; females, 15-49
             )

a <- ABC_sequential(method = "Lenormand",
                    model = get_params,
                    prior = priors,
                    nb_simul = 100,
                    summary_stat_target = sa.prev,
                    p_acc_min = 0.01,
                    n_cluster = 16,
                    use_seed = TRUE)

save(a, file = "/gscratch/csde/kpeebles/abc_het_results_pop_20k.rda")
