## Purpose: Identify single simulation that best matches observed prevalence by least-squares fit to observed data in South Africa from 2002-2012 from among 512 simulations
## Date:    25 March 2018
## Author:  Kathryn Peebles

## Set working directory
setwd("R:/Project/EVONET/user_workspaces/sandbox-kathryn/vacc_cea/abc/het")

## Observed prevalence
sa_year_all     <- c(13, 16, 19, 23)
sa_year_15to24  <- 19
sa_year_sex_age <- c(13, 16, 23)

# Prevalence ages 15-49, observed
sa_prev_15to49   <- c(0.156, 0.162, 0.169, 0.188) # 2002, 2005, 2008, 2012
sa_prev_15to24   <- c(0.087) # 2008
sa_prev_f_15to24 <- c(0.120, 0.169, 0.114) # 2002, 2005, 2012
sa_prev_f_15to49 <- c(0.177, 0.202, 0.232) # 2002, 2005, 2012
sa_prev_m_15to24 <- c(0.061, 0.044, 0.029) # 2002, 2005, 2012
sa_prev_m_15to49 <- c(0.128, 0.117, 0.145) # 2002, 2005, 2012

sa_prev <- c(sa_prev_15to49,
             sa_prev_15to24,
             sa_prev_f_15to24,
             sa_prev_f_15to49,
             sa_prev_m_15to24,
             sa_prev_m_15to49)

## popsumm_frequency parameter is 10 for these simulations
popsumm_freq <- 10

## Specify the time points for each age and sex group for which we'll compare simulated prevalence to observed prevalence
model_year_all     <- sa_year_all * 365/popsumm_freq
model_year_15to24  <- sa_year_15to24 * 365/popsumm_freq
model_year_sex_age <- sa_year_sex_age * 365/popsumm_freq

## Extract prevalence from each model simulation
prev_list <- list()

for(i in 1:8) {
  load(paste0(getwd(), "/sim_sets/sim_set_", i, ".RDATA"))
  model <- eval(parse(text = paste0("sim_set_", i)))
  prev_list[[i]] <- t(sapply(model$popsumm, function(x) {
    rbind(c(x$prev_15to49[model_year_all],
            x$prev_15to24[model_year_15to24],
            x$prev_f_15to24[model_year_sex_age],
            x$prev_f_15to49[model_year_sex_age],
            x$prev_m_15to24[model_year_sex_age],
            x$prev_m_15to49[model_year_sex_age]))
  }))
}

## Combine prevalence from all runs into one table
prev <- do.call(rbind, prev_list)

## Identify simulation with closest fit to observed data by least-squares
sum_squares <- sapply(1:nrow(prev), function(x) {
  sum((sa_prev - prev[x, ]) ^ 2)
})

ls <- which.min(sum_squares)

ls %/% 64 ## simulation ls is in 3rd sim_set
ls %% 64  ## Simulation is #12 in sim_set_3

## Use sim_set_3 for plotting below
model <- sim_set_3

## For below plotting code to work correctly (recycling from earlier approach), reassign ls as # 12
ls <- 12

## Assign prevalence vectors to new objects for code brevity
prev_15to49   <- model$popsumm[[ls]]$prev_15to49
prev_15to24   <- model$popsumm[[ls]]$prev_15to24
prev_f_15to49 <- model$popsumm[[ls]]$prev_f_15to49
prev_f_15to24 <- model$popsumm[[ls]]$prev_f_15to24
prev_m_15to49 <- model$popsumm[[ls]]$prev_m_15to49
prev_m_15to24 <- model$popsumm[[ls]]$prev_m_15to24

## Plot
pdf(file = paste0(getwd(), "/plot_best_fit_ls.pdf"), width = 14, height = 14)

par(mfrow = c(2, 2))

sa_year <- c(seq(1,13,1), 16, 19, 23) * (365/model$param[[1]]$popsumm_frequency)

### Plot of prevalence, ages 15-49
plot(1:length(prev_15to49), 1:length(prev_15to49), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, ages 15-49",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence ages 15-49, observed
sa_prev <- c(0.002, 0.005, 0.010, 0.019, 0.031, 0.048, 0.067, 0.088, 0.108, 0.126, 0.141, 0.153, 0.156, 0.162, 0.169, 0.188)
# CIs for all years except 2002
lb <- c(0.002, 0.005, 0.009, 0.017, 0.029, 0.044, 0.063, 0.084, 0.103, 0.121, 0.136, 0.147, 0.149, 0.155, 0.175)
ub <- c(0.003, 0.006, 0.011, 0.020, 0.033, 0.050, 0.071, 0.092, 0.113, 0.131, 0.146, 0.158, 0.177, 0.184, 0.203)

# Prevalence ages 15-49, model output
lines(x = 1:length(prev_15to49), y = prev_15to49, col = "springgreen4", lwd = 2)

# CIs for 1990-2001, 2005-2012
arrows(sa_year[c(1:12,14:16)], lb, sa_year[c(1:12,14:16)], ub, angle = 90, code = 3, length = 0.05)
# CI for 2002 is estimated from a figure, so should be in grey to indicate it is approximate
arrows(sa_year[13], 0.140, sa_year[13], 0.175, angle = 90, code = 3, length = 0.05, col = "darkgray")

points(x = sa_year[1:12], y = sa_prev[1:12], pch = 19) # Years 1990-2001 from UNAIDS estimates based on ANC prevalence
points(x = sa_year[13:16], y = sa_prev[13:16], pch = 17) # Years from 2002-2012 based on nationally representative surveys (HSRC)

legend("topleft", legend = c("", "", "", "", ""), col = "darkgray", pch = c(rep(NA, 5)), pt.cex = 2, lty = c(NA, NA, 1, NA, NA), bty = "n",
       lwd = 1)

legend("topleft", legend = c("Observed prevalence and 95% CI (UNAIDS estimates from ANC)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys; CI approximate)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "black", "black", "springgreen4"),
       pch = c(19, 17, 17, NA), lty = c(rep(1, 2), NA, 1), lwd = c(rep(1, 3), 2), bty = "n")

# --------------------------------------------------------------------------------------
### Plot of prevalence, ages 15-24
plot(1:length(prev_15to24), 1:length(prev_15to24), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, ages 15-24",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence ages 15-24, model output
lines(x = 1:length(prev_15to24), y = prev_15to24, col = "springgreen4", lwd = 2)

# Prevalence ages 15-24, observed (only 2008, as other years have sex-specific data)
points(x = sa_year[15], y = 0.087, pch = 17)
arrows(sa_year[15], 0.072, sa_year[15], 0.104, angle = 90, code = 3, length = 0.05)

legend("topleft", legend = c("Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "springgreen4"),
       pch = c(17, NA), lty = c(1, 1), lwd = c(1, 2), bty = "n")

# --------------------------------------------------------------------------------------
### Plot of prevalence, females ages 15-24
plot(1:length(prev_f_15to24), 1:length(prev_f_15to24), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, females ages 15-24",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence, females ages 15-24, model output
lines(x = 1:length(prev_f_15to24), y = prev_f_15to24, col = "springgreen4", lwd = 2)

# CIs for 2002 and 2012
arrows(sa_year[13], 0.092, sa_year[13], 0.147, angle = 90, code = 3, length = 0.05)
arrows(sa_year[16], 0.098, sa_year[16], 0.132, angle = 90, code = 3, length = 0.05)
# CI for 2005 is approximate, based on a figure with confidence intervals.
arrows(sa_year[14], 0.145, sa_year[14], 0.20, angle = 90, code = 3, length = 0.05, col = "darkgray")
# Prevalence, females ages 15-24, observed
points(x = c(sa_year[13], sa_year[14], sa_year[16]), y = c(0.120, 0.169, 0.114), pch = 17)

legend("topleft", legend = c("", "", "", ""), col = "darkgray", pch = c(rep(NA, 4)), lty = c(NA, 1, NA, NA), bty = "n", lwd = 1)

legend("topleft", legend = c("Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys; CI approximate)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "black", "springgreen4", "dimgray"),
       pch = c(17, 17, NA), lty = c(1, NA, 1), lwd = c(1, NA, 2), bty = "n")

# --------------------------------------------------------------------------------------
### Plot of prevalence, females ages 15-49
plot(1:length(prev_f_15to49), 1:length(prev_f_15to49), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, females ages 15-49",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence, females ages 15-49, model output
lines(x = 1:length(prev_f_15to49), y = prev_f_15to49, col = "springgreen4", lwd = 2)

# CIs for 2005 and 2012
arrows(sa_year[14], 0.183, sa_year[14], 0.222, angle = 90, code = 3, length = 0.05)
arrows(sa_year[16], 0.213, sa_year[16], 0.251, angle = 90, code = 3, length = 0.05)
# CI for 2002 is approximate, based on a figure with CIs
arrows(sa_year[13], 0.155, sa_year[13], 0.205, angle = 90, code = 3, length = 0.05, col = "darkgray")

# Prevalence, females ages 15-49, observed
points(x = c(sa_year[13], sa_year[14], sa_year[16]), y = c(0.177, 0.202, 0.232), pch = 17)

legend("topleft", legend = c("", "", "", ""), col = "darkgray", pch = c(rep(NA, 4)), lty = c(NA, 1, NA, NA), bty = "n",
       lwd = 1)

legend("topleft", legend = c("Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys; CI approximate)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "black", "springgreen4", "dimgray"),
       pch = c(17, 17, NA), lty = c(1, NA, 1), lwd = c(1, NA, 2), bty = "n")

# --------------------------------------------------------------------------------------
### Plot of prevalence, males ages 15-24
plot(1:length(prev_m_15to24), 1:length(prev_m_15to24), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, males ages 15-24",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence, males ages 15-24, model output
lines(x = 1:length(prev_m_15to24), y = prev_m_15to24, col = "springgreen4", lwd = 2)

# CIs for 2002 and 2012
arrows(sa_year[13], 0.039, sa_year[13], 0.083, angle = 90, code = 3, length = 0.05)
arrows(sa_year[16], 0.021, sa_year[16], 0.039, angle = 90, code = 3, length = 0.05)
# CI for 2005 is approximate, based on a figure
arrows(sa_year[14], 0.030, sa_year[14], 0.065, angle = 90, code = 3, length = 0.05, col = "darkgray")
# Prevalence, males ages 15-24, observed
points(x = c(sa_year[13], sa_year[14], sa_year[16]), y = c(0.061, 0.044, 0.029), pch = 17)

legend("topleft", legend = c("", "", "", ""), col = "darkgray", pch = c(rep(NA, 4)), lty = c(NA, 1, NA, NA), bty = "n", lwd = 1)

legend("topleft", legend = c("Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys; CI approximate)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "black", "springgreen4", "dimgray"),
       pch = c(17, 17, NA), lty = c(1, NA, 1), lwd = c(1, NA, 2), bty = "n")

# --------------------------------------------------------------------------------------
### Plot of prevalence, males ages 15-49
plot(1:length(prev_m_15to49), 1:length(prev_m_15to49), col = "white", xlab = "Year", ylab = "Prevalence", main = "Prevalence, males ages 15-49",
     ylim = c(0.0, 0.35), xaxt = "n")
x_axis <- c(seq(1990,2002,1), 2005, 2008, 2012)
axis(1, at = sa_year, labels = x_axis)

# Prevalence, males ages 15-49, model output
lines(x = 1:length(prev_m_15to49), y = prev_m_15to49, col = "springgreen4", lwd = 2)

# CI for 2002 is approximate, based on figure
arrows(sa_year[13], 0.110, sa_year[13], 0.155, angle = 90, code = 3, length = 0.05, col = "darkgray")
# CIs for 2005 and 2012
arrows(sa_year[14], 0.100, sa_year[14], 0.136, angle = 90, code = 3, length = 0.05)
arrows(sa_year[16], 0.128, sa_year[16], 0.163, angle = 90, code = 3, length = 0.05)

# Prevalence, males ages 15-49, observed
points(x = c(sa_year[13], sa_year[14], sa_year[16]), y = c(0.128, 0.117, 0.145), pch = 17)

legend("topleft", legend = c("", "", "", ""), col = "darkgray", pch = c(rep(NA, 4)), lty = c(NA, 1, NA, NA), bty = "n", lwd = 1)

legend("topleft", legend = c("Observed prevalence and 95% CI (HSRC nationally rep. surveys)",
                             "Observed prevalence and 95% CI (HSRC nationally rep. surveys; CI approximate)",
                             "Model output prevalence of best fit simulation"),
       col = c("black", "black", "springgreen4"),
       pch = c(17, 17, NA), lty = c(1, NA, 1), lwd = c(1, NA, 2), bty = "n")

dev.off()

## Compare prevalence values in dat_saved_12 to confirm that it is the correct saved dat object to use in simulations
load(paste0(getwd(), "/dat_saved_12.RDATA"))

all(dat$popsumm$prev_15to24 == prev_15to24)
all(dat$popsumm$prev_15to49 == prev_15to49)
all(dat$popsumm$prev_f_15to24 == prev_f_15to24)
all(dat$popsumm$prev_f_15to49 == prev_f_15to49)
all(dat$popsumm$prev_m_15to24 == prev_m_15to24)
all(dat$popsumm$prev_m_15to49 == prev_m_15to49)
