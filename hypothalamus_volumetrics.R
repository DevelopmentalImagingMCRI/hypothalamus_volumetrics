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
  
 # Run separately for each region; Superior Tubular provided as example
  
 # Left and Right Tubular Superior 
  
 fit.lts <- brm(bf(supTub_L ~ illdur + sex + age + eTIV, quantile = 0.5), data = hypo_caseonly, family = asym_laplace())
 describe_posterior(fit.lti) %>% print_html()
 plot(conditional_effects(fit.lti))
  
 fit.rts <- brm(bf(supTub_R ~ illdur + sex + age + eTIV, quantile = 0.5), data = hypo_caseonly, family = asym_laplace())
 describe_posterior(fit.rti) %>% print_html()
 plot(conditional_effects(fit.rti))  
  
  
 # ----- PLOTS -----
  
  
 # --- Fatigue Interactions ---
  
 # Plot corrected fatigue interactions for both groups; example using asHyp_R as per paper
 # Replace "asHyp_R" for region/plot of interest
  
 bf.asHyp_R <- stan_glm(asHyp_R ~ group*total_fatigue + sex + age + eTIV,  chains = 10, iter = 5000, warmup = 1000, data = hypo_vols) 
  
 marginal_means <- estimate_means(bf.asHyp_R,at = c("group","total_fatigue")) 
  
 marginal_means %>% ggplot(aes(x = total_fatigue, y= asHyp_R,colour = group)) +
   geom_ribbon(data = marginal_restricted, aes(y = Mean, ymin = CI_low, ymax = CI_high, fill=group), alpha = 0.1, size=0) +
   geom_line(data = marginal_restricted, aes(y = Mean, ymin = CI_low, ymax = CI_high), size = 1.5) +
   xlab("Total Fatigue") +
   ylab("Right Anterior-Superior Volume") +
   theme_minimal() +
   theme(text = element_text(size = 15))
  
  
 # --- Illness Duration ---
  
 # Plot relationship between variables; replace "supTub_R" with ROI 
  
 hypo_caseonly %>% ggplot(aes(x=illdur, y= supTub_R)) +
   geom_point(aes(colour=sex), size = 2) +
   geom_smooth(method ="lm", colour="#666666") +
   ylab("Right Superior Tubular Volume") +
   xlab("Illness Duration") +
   theme_minimal() +
   theme(text = element_text(size = 15)) +
   theme(legend.position="none")
