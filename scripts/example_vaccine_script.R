options(error=traceback)
library(evonet)
###################################

#initial values needed to calculate other parameter values
initial_pop       = 500
years_to_max_coverage = 1#how long to all eligible are vaccinated
max_perc_vaccinated =0.99 #maximum proportion of population eligible for vaccination (.99 upper limit)

#--------------------------------------------------------------
param_list=list(
  nsims            = 1,
  ncores           = 1,
  popsumm_frequency  = 30,
  fast_edgelist      = TRUE,
  save_partner_list  = FALSE,
  #vl_peak_agent_flag = TRUE,
  trans_RR_age       = 1.0,
  min_spvl_allowed = .5,
  n_steps           = 365*30,
  initial_pop       = initial_pop,
  initial_infected  = initial_pop*.10,
  target_stats         = initial_pop*0.7/2,
  #nw_form_terms = "~edges + offset(nodematch('role', diff=TRUE, keep=1:2))"  ,
  vl_peak_agent_flag   = TRUE,  #default FALSE
#Vaccine parameters ----------------------------------------- #
  start_vacc_campaign =  (10*365):(50*365),
  max_perc_vaccinated = max_perc_vaccinated , #maximum percent/proportion of population to be vaccinated
  perc_vaccinated_rate = (max_perc_vaccinated /(1-max_perc_vaccinated ))/(years_to_max_coverage*365),
  vacc_eff_duration = 365*3,
  vaccine_model_id=1, #for Josh and Paul's new model, which one to use, currently only (1) but in future 1,2,3,...
  vacc_trans_prob_decrease =0.8 #proportion (percentage) decrease in trans probs due to vaccine for
                                # vaccine model "1" (baseline vaccine model)
)


evoparams <- do.call(evonet_setup,param_list)



####################################

#add whatever check you want and put in "module" list below
#Example
#test_fxn<-function(dat,at){
  #if(at> (dat$param$start_vacc_campaign+365)){browser()}
#  return(dat)
#}
########################

modules <- c(
   "aging",
  "testing",
  "phi1",
  "update_mu1",
  "covariates1",
  "treatment",
  "viral_update",
  "coital_acts",
  "transmission",
  "evo_departures",
  "evo_arrivals",
  "summary_module")


#estimate network
nw <- nw_setup(evoparams) # Sets up the initial network

#run model
evomodel <- evorun(modules,evoparams,nw)

#aa=traceback()
#aa[1]

#assign model names
model_name = paste("new_vaccine_test1.RData",sep="")
#save model
save(evomodel,file=model_name)
#plot results, saves pdf to working directory, getwd() to see
evoplot(evomodel,name=model_name) #plot to screen (if single sim/core) and write to pdf

##########################


