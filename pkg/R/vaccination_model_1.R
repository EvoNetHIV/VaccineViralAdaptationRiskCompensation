#file contains: phi1
#              update_mu1
#              covariates1
#              marks1
#              draw_m1
#              caclulate_theta1
             
#################

#' @export
phi1 <- function(dat, at) {
  #1) Once vaccination campaign starts, fxn changes vaccination status (attribute "phi") from 
  # NA to "1" to eligibe/selected agents.
  #2) After vacc. campaign has started, will change status of vaccinated agent's whose 
  #efficacy has waned to "0" (stochastic draw)
  #3)this function is put in the "module" list in run script (called each timestep)
  
  
  
  #skip vaccination routine if vacc. campaign hasn't started
  if(at < dat$param$start_vacc_campaign[1]) {return(dat)}
  
  # off/on for already vaccinated
  if(at > dat$param$start_vacc_campaign[1] ) {
    vacc_ix <- which(dat$pop$phi == 1)
    dat$pop$phi[vacc_ix] <- rbinom(length(vacc_ix), 1, 1 - (1/dat$param$vacc_eff_duration))
  }
  
  
  #if vaccination only occurs at discrete intervals, skip rest of routine if vaccination
  #campaign has started but does not occur in this time step
  if(!is.element(at,dat$param$start_vacc_campaign)) {return(dat)}
  #------------------------------------------------
  

  #if designated vacc. level reached (percent of pop vaccianted), don't vacc anymore
  if(length(which(dat$pop$phi == 1 ))/length(which(dat$pop$Status>=0)) > dat$param$max_perc_vaccinated){return(dat)}
  
  
  # Eligible_patients: eligible for care, not vaccinated, not infected
  eligible_index <- which(dat$pop$Status == 0 & 
                            (dat$pop$phi == 0 | is.na(dat$pop$phi)) &
                            dat$pop$eligible_care == 1) 
  
  if(length(eligible_index) == 0) {return(dat)}  #if no agents are eligible
  
  no_vaccinated <- sum(rbinom(length(which(dat$pop$Status>=0)), 1, dat$param$perc_vaccinated_rate)) #denominator is total population alive 
  if(no_vaccinated == 0) {return(dat)}
  
  
  if(no_vaccinated <length(eligible_index)){
    vaccinated_index <- sample(eligible_index, no_vaccinated)
  }else{
    vaccinated_index <- eligible_index
    #if the %coverage in total population alive exceeds #eligible, vaccinate all eligible
  }
  
  dat$pop$phi[vaccinated_index] <- 1
  dat$pop$vacc_init_time[vaccinated_index] <- at
  
  return(dat)
}
#######################

#' @export
update_mu1 <- function (dat,at){
  #code to initialize/update  mu,sigma for each agent or initialize for new agents to model
  #this function is put in the "module" list in run script (called each timestep)
  
  # for model 1, this fxn just returns dat object without modification (no mu,sigma to update)
  
  #note: "m" and "sigma" are agent attributes (see file/fxn "input_parameters_agent_attributes)
  #which means that the default value of the attribute is "NA" until modified
  
  return(dat)
}
#######################

#' @export
covariates1 <- function(dat,at){
  #initializ and update baseline covariates necessary for vaccine model's "calculate_theta" fxn
  #that are not already part of baseline evonet
  #this function is put in the "module" list in run script (called each timestep)
  
  #for model 1, no covariates so just return "dat$covariates" object unmodificed
  
  if(at==2){
    dat$covariates <- NA
  }
  
  return(dat)
}
######################

#' @export
marks1 <- function(dat,index,at){
#for initially infected (model start, at=1), binomial draw to determine whether virus is "sensitive" to vaccine or not
# called in "initialize_infecteds_vl", line 129

#for secondary infections, called in "transmission_bookkeeping_module" , line 206
#for model 1, simply returns marks of infector

  #note at=1 is model initialization
if(at==1){
  rbinom(length(index),1,
         prob=dat$param$perc_virus_vaccine_sens)
  #note: param$perc_virus_vaccine_sens is a user-specified param (default: 0.99)
}else{    
  dat$pop$m[index]
}
  
}
######################
#' @export
draw_m1 <- function(dat){
  # calculate marks, based on  most current values of mu and sigma 
  # called in "calculate_theta1", see below
  
  #for model 1, marks don't change, it is whether infector's virus is sensitive
  #to infectee's vaccine (0/1)
  
  m <- dat$pop$m[dat$infector_id]
  return(m)
}
#######################

#' @export
calculate_theta1 <- function(dat){
  # draw updated marks (of virus)
  # then caclulate theta based on dat$m, dat$phi, dat$covariates
  # output theta to adjust raw trans probs
  #called at "transmission_main_module" line 157
  
  #for model 1, mark is whether virus of infector is sensitive to vaccine of
  # infectee

  m <- draw_m1(dat)
  #theta is  vaccination status (0/1) * 
  #             virus susceptibility to vaccine (0/1) *
  #                percent reduction due to vaccine in trans. prob.
  
  theta <- rep(1,length(dat$susceptible_id))
  index <- which(dat$pop$phi[dat$susceptible_id]==1)
  if(length(index)==1){
    theta[index]  <-(dat$pop$phi[dat$susceptible_id])[index]*m[index]*dat$param$vacc_trans_prob_decrease
  }
  return(theta)
}


