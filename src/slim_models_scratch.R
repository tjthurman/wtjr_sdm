library(tidyverse)


change.df %>% 
  filter(Long > -114) %>% 
  slice_max(change_probBrown, n = 200)

rpois(10, 4)
  

current.df <- as(current_pheno, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1]))




# Rootsolve for a

fitness_at_mismatch <- function(sd) {
  
  return(dnorm(1, 0, sd)/dnorm(0, 0, sd))
  
}

fitness_at_mismatch(0.6)

diff_from_target  <- function(sd, target, annual_survival) {
  
  
  return(fitness_at_mismatch(sd) - target)
}

uniroot(f = diff_from_target, lower = 0.1, upper = 1, target = 0.520)


# Kardos logistic growth

kardos_log_growth <- function(mean_fitness, K, N) {
  
  return(exp(log(mean_fitness)*(1- (N/K))))
}

my_log_growth <- function(mean_fitness, K, N) {
  
  mean_fitness * K/N

}


kardos_log_growth(0.9, 5000, 2500)
my_log_growth(0.9, 5000, 2500)

x <- seq(-5, 5, 0.01)

multiplied <- dnorm(x = x, mean = 0, 1)/dnorm(0, mean = 0, 1)
added <- dnorm(x = x, mean = 0, 1) + (1 - dnorm(0, 0,1))

sd(multiplied)
sd(added)


plot(x = x, y = multiplied, type = "l", ylim = c(0,1))
lines(x = x, y = added, col = "red")
lines(x = x, y = dnorm(x), col = "blue")
lines(x = x, y = gauss_vec(x, 1, 0, 1), col = "orange")
abline(v = -1, lty = "dotted")


(dnorm(x = 1, mean = 0, 1)/dnorm(0, mean = 0, 1))
(dnorm(x = 1, mean = 0, 1) + (1 - dnorm(0, 0,1)))

gaussian_func <- function(x,a,b,c) {
  
  a * exp(-1 *((x-b)^2)/(2*c^2))
}

gaussian_func(1, 1, 0, 0.65)

gauss_vec <- Vectorize(gaussian_func)



diff_from_target  <- function(c, a, b, x, target) {
  
  gauss_vec(x = x, a = a, b = b, c = c) - target
}

# Find sd for a given # of weeks of winter
uniroot(f = diff_from_target, lower = 0, upper = 1, a = 1, b = 0, x = -1, target = 0.85^16)


survival_rate = 0.2
K = 50000

to_replace <- function(mean_children, min_surv, K) {

  (min_surv*K)/2*mean_children + (min_surv*K) - K
      
}


x <- seq(0.05, 0.5, by = 0.025)
y <- NULL
for (z in x) {
  y <- c(y, uniroot(f = to_replace, lower = 0, upper = 1000, min_surv = z, K = 1000)$root)
}

plot(x, y)



tibble(x, y)

2^(0.5)

uniroot(f = to_replace, lower = 0, upper = 50, min_surv = 0.2, K = 1000)$root

new_pop_size(0.2, 40000, 6)


rpois(100, 4.2)
