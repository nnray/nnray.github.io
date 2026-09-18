# simulation function ##########################################################

simulation <- function(n, corr, beta_1, sims){
  
  library(MASS)
  
  set.seed(1)

  XZ <- mvrnorm(n = n,
                   mu = c(0, 0),
                   Sigma = matrix(c(1, corr, corr, 1), ncol = 2),
                   empirical = T)
  
  X = XZ[, 1]
  Z = XZ[, 2]
  
  MZ = diag(n) - (Z %*% solve(t(Z) %*% Z) %*% t(Z))
  
  results <- data.frame(
    sample_size = n,
    b1 = rep(NA, n),
    R_squared = rep(NA, n),
    VIF = rep(NA, n),
    p_value = rep(NA, n)
  )
  
  for(i in 1:sims){
    
    results[i, "sample_size"] <- n
    
    Y = beta_1 * X + 1 * Z + rnorm(n)
    
    results[i, "b1"] <- solve(t(X) %*% MZ %*% X) %*% (t(X) %*% MZ %*% Y)
    
    R_sq = 1 - ((t(X) %*% MZ %*% X) / (t(X) %*% X))
    
    results[i, "R_squared"] <- R_sq
    
    results[i, "VIF"] <- 1 / (1 - R_sq)
    
    results[i, "p_value"] <- summary(lm(Y ~ -1 + X + Z))[["coefficients"]]["X", "Pr(>|t|)"]
    
  }
  
  return(results)
  
}

# executing simulation function ################################################

out_0.99 <- rbind(
  simulation(n = 50, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 100, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 150, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 200, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 250, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 300, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 350, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 400, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 450, corr = 0.99, beta_1 = 1, sims = 1000),
  simulation(n = 500, corr = 0.99, beta_1 = 1, sims = 1000)
)

out_0.90 <- rbind(
  simulation(n = 50, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 100, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 150, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 200, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 250, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 300, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 350, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 400, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 450, corr = 0.9, beta_1 = 1, sims = 1000),
  simulation(n = 500, corr = 0.9, beta_1 = 1, sims = 1000)
)

# checking VIFs
out_0.99$VIF |> unique()
out_0.90$VIF |> unique()

# plotting output ##############################################################

library(tidyverse)

out_0.99 <- out_0.99 |> 
  mutate(reject = ifelse(p_value <= 0.05, 1, 0)) |>
  group_by(sample_size) |> 
  mutate(rate = mean(reject)) |> 
  ungroup()

out_0.90 <- out_0.90 |> 
  mutate(reject = ifelse(p_value <= 0.05, 1, 0)) |>
  group_by(sample_size) |> 
  mutate(rate = mean(reject)) |> 
  ungroup()

pdf("corr_0.99.pdf", width = 5, height = 5)

ggplot(out_0.99) +
  geom_point(aes(x = sample_size, y = rate)) +
  theme_bw() +
  labs(x = "Sample Size", y = "Power",
       title = bquote("Correlation between" ~ italic("X") ~ "and" ~ italic("Z") ~ "equals 0.99")) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

dev.off()

pdf("corr_0.90.pdf", width = 5, height = 5)

ggplot(out_0.90) +
  geom_point(aes(x = sample_size, y = rate)) +
  theme_bw() +
  labs(x = "Sample Size", y = "Power",
       title = bquote("Correlation between" ~ italic("X") ~ "and" ~ italic("Z") ~ "equals 0.90")) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

dev.off()
