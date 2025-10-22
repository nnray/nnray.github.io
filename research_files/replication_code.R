# defining the experiment ######################################################

simulation <- function(
    
  units, periods, # N and T
  y_autocorrelation, # autocorrelation in the outcome
  treatment_effect, treatment_period, # true short-run effect and when it occurs
  sims, return_data, # number of simulations and option to return the data instead
  burn_periods # number of time points to burn
  
){
  
  # required libraries
  library(tidyverse); library(fixest); library(did); library(foreach);
  library(doRNG); library(doParallel)
  
  set.seed(1)
  
  data <- data.frame(
    n = rep(1:units, each = (burn_periods + periods)),
    t = rep(1:(burn_periods + periods), units),
    y = NA
  )
  
  # half of the units are randomly selected to receive treatment
  treated_units <- sample(units, 0.5 * units)
  data <- mutate(data, treat = ifelse(n %in% treated_units, 1, 0))
  
  # creating a variable for the true treatment for easier simulating
  data <- mutate(data, treatment_effect = ifelse(treat == 1 & t == (burn_periods + treatment_period), treatment_effect, 0))
  
  # creating the "post" variable that will be used in the DiD regression (using `fixest`)
  data <- mutate(data, d = ifelse(treat == 1 & t >= (burn_periods + treatment_period), 1, 0))
  
  # creating the "first treat" variable used in `did` package estimation
  data <- mutate(data, first.treat = ifelse(treat == 1, (burn_periods + treatment_period), 0))
  
  # empty data.frame to store results
  results <- data.frame(
    
    # number of units
    units = NA,
    
    # estimates of the "post" variable d
    d_hat_iid = NA, d_hat_clustered = NA, d_hat_ldv = NA, d_hat_adl11 = NA,
    
    # estimates of the standard deviation of the "post" variable d
    se_hat_iid = NA, se_hat_clustered = NA, se_hat_ldv = NA, se_hat_adl11 = NA,
    
    # p-values for the "post" variable d
    p_iid = NA, p_clustered = NA, p_ldv = NA, p_adl11 = NA
    
  )
  
  # number of cores to use for parallelization and registering them
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  # parallelized simulation
  results <- foreach(
    s = 1:sims, .combine = rbind, .packages = c("tidyverse", "fixest")
    ) %dorng% {
    
    for(i in 1:units){
      
      # within each unit i, making the first observation a random draw (iid across units)
      data[, "y"][data[ , "n"] == i][1] <- rnorm(1)
      
      # making y autocorrelated within units, plus treatment (non iid within units)
      for(j in 2:(burn_periods + periods)){
        
        data[, "y"][data[ , "n"] == i][j] <- y_autocorrelation *
          data[, "y"][data[ , "n"] == i][j - 1] + rnorm(1) +
          data[, "treatment_effect"][data[ , "n"] == i][j]
        
      }
      
    }
    
    data <- filter(data, t > burn_periods)
    
    results[s, "units"] <- units
    
    # calculating lags of y and d for ldv and adl models
    lag_data <- data %>% group_by(n) %>% mutate(y_lag = lag(y), d_lag = lag(d))
    
    # static model with standard errors estimated iid, TWFE
    iid <- feols(y ~ d | n + t, data, vcov = "iid")
    results[s, "d_hat_iid"] <- iid[["coefficients"]][["d"]]
    results[s, "se_hat_iid"] <- iid[["se"]][["d"]]
    results[s, "p_iid"] <- iid[["coeftable"]][4]
    
    # static model with standard errors clustered on unit, TWFE
    clustered <- feols(y ~ d | n + t, data, vcov = cluster ~ n)
    results[s, "d_hat_clustered"] <- clustered[["coefficients"]][["d"]]
    results[s, "se_hat_clustered"] <- summary(clustered)[["se"]][["d"]]
    results[s, "p_clustered"] <- summary(clustered)[["coeftable"]][4]
    
    # ldv model with standard errors estimated iid, TWFE
    ldv <- feols(y ~ y_lag + d | n + t, lag_data, vcov = "iid")
    results[s, "d_hat_ldv"] <- ldv[["coefficients"]][["d"]]
    results[s, "se_hat_ldv"] <- ldv[["se"]][["d"]]
    results[s, "p_ldv"] <- ldv[["coeftable"]][2, 4]
    
    # adl(1,1) model with standard errors estimated iid, TWFE
    adl11 <- feols(y ~ y_lag + d + d_lag | n + t, lag_data, vcov = "iid")
    results[s, "d_hat_adl11"] <- adl11[["coefficients"]][["d"]]
    results[s, "se_hat_adl11"] <- adl11[["se"]][["d"]]
    results[s, "p_adl11"] <- adl11[["coeftable"]][2, 4]
    
    results[s, ]
    
  }
  
  # close back-end
  stopCluster(cl)
  
  if(return_data == T){return(data)} else{return(results)}
  
}

# conducting the experiments ###################################################

# simulation 1 (used for Figure 1)
system.time(output <- simulation(
  units = 50,
  periods = 10,
  y_autocorrelation = 1,
  treatment_effect = 0,
  treatment_period = 5,
  sims = 1000,
  return_data = F,
  burn_periods = 100
))
# took 179.82 seconds (2.997 minutes)

# simulation 2 (used for Figure 2)
# system.time(output <- simulation(
#   units = 50,
#   periods = 10,
#   y_autocorrelation = 1.05,
#   treatment_effect = 0,
#   treatment_period = 5,
#   sims = 1000,
#   return_data = F,
#   burn_periods = 100
# ))
# took a comparable amount of time relative to simulation 1

# simulation 3 (used for Figure 3)
# system.time(output <- simulation(
#   units = 1000,
#   periods = 10,
#   y_autocorrelation = 1,
#   treatment_effect = 0,
#   treatment_period = 5,
#   sims = 1000,
#   return_data = F,
#   burn_periods = 100
# ))
# took 7056.27 seconds, or 117.6 minutes

# calculating and plotting results #############################################

# plotting estimates of the "post" variable d
estimates <- ggplot(output) +
  geom_point(aes(x = 0, y = d_hat_iid)) +
  geom_point(aes(x = 1, y = d_hat_clustered)) +
  geom_point(aes(x = 2, y = d_hat_ldv)) +
  geom_point(aes(x = 3, y = d_hat_adl11)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2, 3),
                     labels = c("Static; iid se", "Static; clustered se", "LDV; iid se", "ADL(1,1); iid se")) +
  labs(title = "Estimated Treatment Effects (True Effect is 0)", x = "Estimators", y = "Effect") +
  theme_bw()

# calculating the standard deviation of d_hat
output <- mutate(output,
                 sd_clustered = sd(d_hat_clustered),
                 sd_iid = sd(d_hat_iid),
                 sd_ldv = sd(d_hat_ldv),
                 sd_adl11 = sd(d_hat_adl11))

# plotting estimated standard deviations of d_hat
ses <- ggplot(output) +
  geom_point(aes(x = 0, y = se_hat_iid)) +
  geom_point(aes(x = 0, y = sd_iid, colour = "SD of Effect Estimates")) +
  geom_point(aes(x = 1, y = se_hat_clustered)) +
  geom_point(aes(x = 1, y = sd_clustered, colour = "SD of Effect Estimates")) +
  geom_point(aes(x = 2, y = se_hat_ldv)) +
  geom_point(aes(x = 2, y = sd_ldv, colour = "SD of Effect Estimates")) +
  geom_point(aes(x = 3, y = se_hat_adl11)) +
  geom_point(aes(x = 3, y = sd_adl11, colour = "SD of Effect Estimates")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3),
                     labels = c("Static; iid se", "Static; clustered se", "LDV; iid se", "ADL(1,1); iid se")) +
  labs(title = "Estimated Standard Deviations of Effects", colour = "Legend", x = "Estimators", y = "SE") +
  theme_bw() +
  theme(legend.position = "bottom")

# calculating false positive rates using p-values
output <- mutate(output,
                 false_positive_iid = ifelse(p_iid <= 0.05, 1, 0),
                 rate_iid = sum(false_positive_iid) / nrow(output),
                 false_positive_clustered = ifelse(p_clustered <= 0.05, 1, 0),
                 rate_clustered = sum(false_positive_clustered) / nrow(output),
                 false_positive_ldv = ifelse(p_ldv <= 0.05, 1, 0),
                 rate_ldv = sum(false_positive_ldv) / nrow(output),
                 false_positive_adl11 = ifelse(p_adl11 <= 0.05, 1, 0),
                 rate_adl11 = sum(false_positive_adl11) / nrow(output))

# plotting false positive rates
rates <- ggplot(output) +
  geom_point(aes(x = 0, y = rate_iid)) +
  geom_point(aes(x = 1, y = rate_clustered)) +
  geom_point(aes(x = 2, y = rate_ldv)) +
  geom_point(aes(x = 3, y = rate_adl11)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2, 3),
                     labels = c("Static; iid se", "Static; clustered se", "LDV; iid se", "ADL(1,1); iid se")) +
  labs(title = "False Positive Rates", x = "Estimators", y = "False Positive Rate") +
  theme_bw()

# combining and saving plots ###################################################

library(cowplot); library(here)

# results when y  is a unit root (simulation 1)
pdf(here("Unit_Root.pdf"), width = 15, height = 8.75)

plot_grid(
  estimates, ses, rates,
  nrow = 1
)

dev.off()

# results when y is explosive (simulation 2)
# pdf(here("Explosive.pdf"), width = 15, height = 8.75)
# 
# plot_grid(
#   estimates, ses, rates,
#   nrow = 1
# )
# 
# dev.off()

# results for unit root and 1000 units (simulation 3)
# pdf(here("Unit_Root_1000_units.pdf"), width = 15, height = 8.75)
# 
# plot_grid(
#   estimates, ses, rates,
#   nrow = 1
# )
# 
# dev.off()
