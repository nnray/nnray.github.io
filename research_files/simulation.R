library(MASS)

set.seed(1)

simulation <- function(n, corr, beta_1, sims){

  X1_X2 <- mvrnorm(n = n,
                   mu = c(0, 0),
                   Sigma = matrix(c(1, corr, corr, 1), ncol = 2),
                   empirical = T)
  
  X1 = X1_X2[, 1]
  X2 = X1_X2[, 2]
  
  results <- data.frame(
    sample_size = n,
    estimate = rep(NA, n),
    R_squared = rep(NA, n),
    VIF = rep(NA, n),
    p_value = rep(NA, n)
  )
  
  for(i in 1:sims){
    
    results[i, "sample_size"] <- n
    
    Y = beta_1 * X1 + 1 * X2 + rnorm(n)
    
    M2 = diag(n) - (X2 %*% solve(t(X2) %*% X2) %*% t(X2))
    
    results[i, "estimate"] <- solve(t(X1) %*% M2 %*% X1) %*% (t(X1) %*% M2 %*% Y)
    
    R_sq = 1 - ((t(X1) %*% M2 %*% X1) / (t(X1) %*% X1))
    
    results[i, "R_squared"] <- R_sq
    
    results[i, "VIF"] <- 1 / (1 - R_sq)
    
    results[i, "p_value"] <- summary(lm(Y ~ -1 + X1 + X2))[["coefficients"]]["X1", "Pr(>|t|)"]
    
  }
  
  return(results)
  
}

out <- rbind(
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

out$VIF |> unique()

library(tidyverse)

ggplot(out) +
  geom_point(aes(x = sample_size, y = p_value)) +
  theme_bw()

out <- out |> 
  mutate(reject = ifelse(p_value <= 0.05, 1, 0)) |>
  group_by(sample_size) |> 
  mutate(rate = mean(reject))

ggplot(out) +
  geom_point(aes(x = sample_size, y = rate)) +
  theme_bw() +
  labs(x = "Sample Size", y = "Power",
       title = bquote("Correlation between" ~ bold("X") ~ "and" ~ bold("Z") ~ "equals 0.9"))
