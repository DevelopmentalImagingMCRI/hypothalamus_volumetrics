#
#    CODE FOR HYPOTHALAMUS VOLUMETRIC ANALYSIS 
#
#  Created by Hollie Byrne & Marc Seal, October 2022
#     contact : hollie@student.unimelb.edu.au
#
 
 
# code in chunks to call from R Markdown document
 
 
# ---- LoadPackages ----
 
library(tidyverse) 
library(easystats) 
library(broom)         # marginal effects and means
library(broom.mixed)   # as above
library(emmeans)       # as above
library(gtsummary)     # for html tables
library(ggeffects)     # for calculating and plotting marginal means 
library(rstanarm)      # for Bayes Analysis 
library(BayesFactor)   # as above
library(brms)          # as above
 
 
## ---- LoadData ---- 
 
# customise path to data
  hypo_vols <- read_csv(file="$path_to_data")

# add labels to categorical variables
  hypo_vols$group <- fct_recode(as.factor(hypo_vols$group), control = "1", case = "0") 
  hypo_vols$sex <- fct_recode(as.factor(hypo_vols$sex), female = "1", male = "0")

# view df and check variable structure 
  glimpse(hypo_vols) 
 
 
# ---- CreateLongdf ----

# create a "nested long Dataframe" by region volumes

df.long <- hypo_vols %>%
          pivot_longer(cols = c(aiHyp_L, asHyp_L, posHyp_L,	infTub_L,	supTub_L,	aiHyp_R, asHyp_R,	posHyp_R,	infTub_R,	supTub_R,	whole_L,	whole_R),
                        names_to = "region", values_to = "score")
 
df.long.nest <- nest(group_by(df.long, region))
 
 
# ---- Aim1-CompareVolumes ----
 
# Bayesian linear model (continuous variable) - can group status predict regional volume?
 
df.long.nest <- mutate(df.long.nest, stan_model = map(data, ~stan_glm(score ~ group + sex + age + eTIV, chains = 10, iter = 5000, 
                                                                      warmup = 1000, data = .x))) 
 
df.long.nest <- df.long.nest %>%
  mutate(stan_model_summary = map(stan_model,summary)) %>% 
  mutate(stan_tidied = map(stan_model, tidy)) %>% 
  mutate(bf.list = map(stan_model, ~ bayesfactor_parameters(.x, null=0)))
 
bf.tidy <- df.long.nest %>% dplyr::select(region,bf.list) %>%
  unnest(cols=bf.list) %>%
  mutate(BF = exp(log_BF))
 
bf.tidy %>% filter(Parameter=="groupcontrol") %>%
  dplyr::select(region,log_BF,BF) %>% print_html()
 
df.long.nest$stan_model # view output for all models 
 
report(df.long.nest$stan_model[[1]]) # the number here responds to which regional hypothalamus volume in the df.long.nest list; 
                                     # run to view for each model (region) of interest
 
describe_posterior(df.long.nest$stan_model[[1]]) # output for first regional hypothalamus volume
 
 
# ---- Aim2-FatigueInteractions ----
 
# Bayesian linear models (continuous variable) - can regional volumes predict fatigue scores by group?  
 
hypo_caseonly <- hypo_vols %>% filter(group =="case") # Run separately for controls by replacing with `filter(group=="control")` and following steps
 
df.long.caseonly <- hypo_caseonly %>%
  pivot_longer(cols = c(aiHyp_L, asHyp_L, posHyp_L,	infTub_L,	supTub_L,	aiHyp_R, asHyp_R,	posHyp_R,	infTub_R,	supTub_R,	whole_L,	whole_R),
               names_to = "region", values_to = "score")
 
df.long.nest.caseonly <- nest(group_by(df.long.caseonly, region))
 
df.long.nest.caseonly <- mutate(df.long.nest.caseonly, stan_model = map(data, ~stan_glm(score ~ total_fatigue + sex + age + eTIV,  chains = 10, iter = 5000, 
                                                                          warmup = 1000, data = .x))) 
  
 df.long.nest.caseonly <- df.long.nest.caseonly %>%
   mutate(stan_model_summary = map(stan_model,summary)) %>% 
   mutate(stan_tidied = map(stan_model, tidy)) %>% 
   mutate(bf.list = map(stan_model, ~ bayesfactor_parameters(.x, null=0))) 
  
 bf.tidy <- df.long.nest.caseonly %>% dplyr::select(region,bf.list) %>%
   unnest(cols=bf.list) %>%
   mutate(BF = exp(log_BF))
  
 bf.table <- bf.tidy %>% 
   dplyr::select(region,BF) %>% 
   print_html() 
  
 print(bf.table) # BF values printed in order of (Intercept), total_fatigue, sex, age and eTIV 
  
 df.long.nest.caseonly$stan_model # view output for all models 
  
 describe_posterior(df.long.nest.caseonly$stan_model[[1]]) # where 1 is replaced by ROI
  
 # Repeat for control group 
  
  
 # ---- Aim3-IllnessDuration ----
  
 # Bayesian quantile regression model (ordinal variable) - is there a relationship between order of illness duration categories and regional hypothalamus volumes?
  
 # Reference: Fit Bayesian Generalized (Non-)Linear Multivariate Multilevel Models using the Asymmetric Laplace Distribution
 # https://paul-buerkner.github.io/brms/reference/brm.html

  
 # create illness as a factor variable   
 hypo_caseonly$ill.factor <- as_factor(hypo_caseonly$illdur)
 class(hypo_caseonly$ill.factor)

 # For ordered factors, create illdurO or illness duration ordered 

 hypo_caseonly <- mutate(hypo_caseonly, illdurO = ordered(illdur))


 #  check class:
 hypo_caseonly$illdurO
 class(hypo_caseonly$illdurO)

 hypo_caseonly <- hypo_caseonly %>%
   mutate(age.scale = scale(age,center = TRUE, scale = TRUE),
         aiHyp_L_scale = scale(aiHyp_L,center = TRUE, scale = TRUE),
         aiHyp_R_scale = scale(aiHyp_R,center = TRUE, scale = TRUE)) # and so on ...

 # run for each region, then run analysis:

 fit.aiHyp_L <- brm(deg_aiHyp_L_scale ~ mo(illdurO) + sex + age.scale, iter= 10000, data = hypo_caseonly)
 describe_posterior(fit.aiHyp_L) %>% print_html()
 plot(conditional_effects(fit.aiHyp_L))

 fit.aiHyp_R <- brm(deg_aiHyp_R_scale ~ mo(illdurO) + sex + age.scale, iter= 10000, data = hypo_caseonly)
 describe_posterior(fit.aiHyp_R) %>% print_html()
 plot(conditional_effects(fit.aiHyp_R)) 
  
