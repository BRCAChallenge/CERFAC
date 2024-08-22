##Version 1.0.
## still need to adjust OR by ethnicity

library(tidyverse)

##This script takes case control observational data as input and computes odds ratios
# for user defined analytic subsets. It also calculates the proportion of variants
# that are pathogenic in each subset.

#Lines 11-102 define functions that are used.
#Lines 104-118 set user defined thresholds
#Lines 121-209 do computations

#functions for proportion pathogenic calculation
probfromrelativerisk = function(delta0, rho)
{
  rho * delta0
}

probfromoddsratio = function(delta0, phi)
{
  phi * delta0 / ( 1 + delta0 * (phi-1))
}

delta = function(p,l,d0,d1)
{
  d0 * exp(-p*l) + d1*(1-exp(-p*l))
}

alpha = function(p,l,d0,d1)
{
  1 - (1-d0)*exp(-l) / (1-delta(p,l,d0,d1))
}

beta = function(p,l,d0,d1)
{
  1 - d0*exp(-l) / delta(p,l,d0,d1)
}

lambda = function(a,b,d0)
{
  log( ((1-b)*(1-d0) + (1-a)*d0) / (1-a)/(1-b))
}

pp = function(a,b,d0,d1)
{
  p = (d1-d0) * ((1-b)*(1-d0) + (1-a)*d0)
  p = p / ( (1-b)*d1*(1-d0) - (1-a)*d0*(1-d1) )
  log(p) / lambda(a,b,d0)
}

propallpath = function(x,n,y,m,d0,d1)
{
  p = pp(x/n,y/m,d0,d1)
  p[p<0] = 0
  p[p>1] = 1
  p
}

propobspath = function(x,n,y,m,d0,d1)
{
  p = propallpath(x,n,y,m,d0,d1)
  l = lambda(x/n,y/m,d0)
  d = delta(p,l,d0,d1)
  
  ai = (p * (1-d1) + (1-p) * (1-d)) / (1-d)
  bi = (p * d1 + (1-p)*d) / d 
  pobs = n * ai + m * bi
  aip = p * (1-d1) /(1-d)
  bip = p * d1 / d
  ppathobs = n*aip + m * bip
  ppathobs/pobs
}

bootall = function(x,n,y,m,d0,d1,b=1000)
{
  xjs = rbinom(b,n,x/n)
  yjs = rbinom(b,m,y/m)
  propallpath(xjs,n,yjs,m,d0,d1)
}

bootobs = function(x,n,y,m,d0,d1,b=1000)
{
  xjs = rbinom(b,n,x/n)
  yjs = rbinom(b,m,y/m)
  propobspath(xjs,n,yjs,m,d0,d1)
}

loglike = function(x,n,y,m,p,l,d0,d1)
{
  a = alpha(p,l,d0,d1)
  b = beta(p,l,d0,d1)
  x * log(a) + (n-x)*log(1-a) + y * log(b) + (m-y) * log(1-b)
}

llplot = function(x,n,y,m,d0,d1)
{
  phat = propallpath(x,n,y,m,d0,d1)
  allp = (0:100)/100
  lhat = lambda(x/n,y/m,d0)
  ll = loglike(x,n,y,m,allp,lhat,d0,d1)
  xl = paste("MLE for p is ",round(100*phat,2),"%",sep="")
  plot(allp,ll,type="l",xlab=xl,ylab="Log likelihood")
  lines(c(phat,phat),range(ll),col=2)
}

#define functional thresholds
loss_of_function_threshold = 0.9
functional_threshold = 0.1

#define total case and control counts
control = 807162
case = 100000

#define population frequency of disease and odds ratio thresholds
d0 = 0.12
oddsrat = 5
d1 = probfromoddsratio(d0,oddsrat)

#define location of file
input_data = '/Users/scottpew/Desktop/PhD/Lab/cerfac/BRCA1_variants_clinical.csv'


#import data
test_data = read.csv(input_data)

#create new column for functional classification
# 0 = non-carrier
# 1 = functional
# 2 = indeterminate
# 3 = loss of function
test_data$func_class = ifelse(test_data$score >= loss_of_function_threshold, 3,
                              ifelse(test_data$score <= functional_threshold, 1,
                                     2))
#sum gnomAD observations
test_data$control_obs = rowSums(test_data[, c(27,31,35,39,43,47,51,55,59,
                                              66,70,74,78,82,86,90,94,98)],
                                na.rm =T)
#keep unique instances of variants
test_data_unique = test_data %>%
  distinct(hgvs_nt, .keep_all = T)

#summarize
data_summary = test_data_unique %>%
  filter(!is.na(func_class)) %>%
  group_by(func_class) %>%
  summarise(control_obs = sum(control_obs),
            case_obs = sum(OBS))

non_carriers = data.frame(
  func_class = 0,
  control_obs = control - sum(data_summary$control_obs),
  case_obs = case - sum(data_summary$case_obs)
)

data_summary_final = rbind(non_carriers, data_summary)

#expand to long table for logistic regression
long_data = data_summary_final %>%
  gather(key = "Category", value = "Count", -func_class) %>%
  rowwise() %>%
  do(data.frame(status = rep(.$Category, .$Count),
                func_class = rep(.$func_class, .$Count)))

# Convert status to a factor
# 0 = controls
# 1 = cases
long_data$status = ifelse(long_data$status == "control_obs", 0,
                          1)
long_data$status = factor(long_data$status)
long_data$func_class = factor(long_data$func_class)

#write.csv(long_data,"/Users/scottpew/Desktop/PhD/Lab/cerfac/test_or.csv")

#glm
model = glm(status ~ factor(func_class), data = long_data, family = binomial)
summary(model)

#extract coefficients and standard errors
coef_estimates = summary(model)$coefficients
log_odds = coef_estimates[, "Estimate"]
se = coef_estimates[, "Std. Error"]

#calculate 95% CI
z_value = qnorm(0.975) #1.96 for 95%CI (finds the critical value in the upper tail (0.05/2 = 0.025. 1-0.25 = 0.975))
ci_log_odds_lower = log_odds - z_value * se
ci_log_odds_upper = log_odds + z_value * se

#convert to odds ratios and corresponding CIs
odds_ratios = exp(log_odds)
ci_odds_ratios_lower = exp(ci_log_odds_lower)
ci_odds_ratios_upper = exp(ci_log_odds_upper)

#make new df
or_results = data.frame(
  OR = odds_ratios,
  CI_lower = ci_odds_ratios_lower,
  CI_upper = ci_odds_ratios_upper
)

or_results = or_results %>%
  mutate(across(where(is.numeric), ~round(., 2)))


prop_path = cbind(data_summary_final, or_results)
prop_path$controls = control
prop_path$cases = case
prop_path$pop_freq = d0
prop_path$prop_threshold = d1
prop_path$prop_path = apply(prop_path, 1, function(row) propallpath(row[2], row[7], row[3], row[8], row[9], row[10]))
prop_path$prop_path = round(prop_path$prop_path, 2)
print(prop_path)
