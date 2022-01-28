# Code for figuring out the standard deviation of
# the gaussian fitness function that corresponds to various
# weekly reductions in survival probability. 

# First, a simple Gaussian function:
gaussian_func <- function(x,a,b,c) {
  
  a * exp(-1 *((x-b)^2)/(2*c^2))
}
gauss_vec <- Vectorize(gaussian_func)

# In our case, we want:
# a = 1 (probability of survival is 100% if perfectly matched) 
# b = 0 (maximize the gaussian function when mismatch = 0)
# And we're trying to figure out what C to use. 

# Next, we can set up a function to use for our root solver:
# This will work by altering the possible c values, finding the one
# that will make our gaussian function equal our target value. 
diff_from_target  <- function(c, a, b, x, target) {
  
  gauss_vec(x = x, a = a, b = b, c = c) - target
}

# Next, we want to decide on target values. 
# We're dealing with this in terms of reduction in weekly survival probabilities,
# to match the empirical estimates for SSH from Zimova et al. 2016
# That paper estimates a 7% weekly reduction in survival for completely mismatched hares. 
# We'll try 7%, as well as 5% and 10% to consider stronger and weaker selection.
# At the locations we're simulating, there are currently 16 weeks of snow cover,
# So there are 16 weeks for mismatch to affect weekly survival. . 
# We will assume survival is guaranteed outside of this time (modelling only mismatch selection).
# We will also assume that perfectly-matched hares have 100% survival during the winter (as we mention above).
# So, a 5% reduction in weekly survival would be a 95% survival rate.
# 7%  reduction in weekly survival = 93% survival rate
# And a 10% reduction in weekly survival = 90% survival rates. 
# To get annual survival rates from these weekly survival rates, we just take them to the 16th power,
# since that's the number of weeks of winter
sd_5perc <- uniroot(f = diff_from_target, lower = 0, upper = 1, a = 1, b = 0, x = -1, target = 0.95^16)$root
sd_7perc <- uniroot(f = diff_from_target, lower = 0, upper = 1, a = 1, b = 0, x = -1, target = 0.93^16)$root
sd_10perc <- uniroot(f = diff_from_target, lower = 0, upper = 1, a = 1, b = 0, x = -1, target = 0.90^16)$root

# We can double check that these lead to the values we want
# First, survival percentages should be 1 with no mismatch
gaussian_func(x = 0, a = 1, b = 0, c = sd_5perc)
gaussian_func(x = 0, a = 1, b = 0, c = sd_7perc)
gaussian_func(x = 0, a = 1, b = 0, c = sd_10perc)
# Then, survival percentages should equal our target for hares that are completely mismatched due to being 
# completely brown when they should be white
# Due to floating issues, these won't match exactly, but should be quite close
all.equal(gaussian_func(x = 1, a = 1, b = 0, c = sd_5perc), 0.95^16)
all.equal(gaussian_func(x = 1, a = 1, b = 0, c = sd_7perc), 0.93^16)
all.equal(gaussian_func(x = 1, a = 1, b = 0, c = sd_10perc), 0.90^16)
# Looks good. For the simulations I round to 4 decimal places anyway, so this slight numerical imprecision doesn't matter. 


# In retrospect, since a = 1 and b = 0, we can skip all the root solving and just solve for c:
find_c <- function(target_surv_rate) {
  sqrt(-1/(2*log(target_surv_rate)))
}

all.equal(find_c(0.95^16), sd_5perc)
all.equal(find_c(0.93^16), sd_7perc)
all.equal(find_c(0.90^16), sd_10perc)
# And you get the same answers