## Summarized script for "Targeted therapy using Polymyxin B hemadsorption in patients with sepsis: A post-hoc analysis of the JSEPTIC-DIC study and the EUPHRATES trial"
## Heterogeneous treatment effects of PMX-HA in patients with sepsis using Jseptic-DIC cohort
## Analyzed by Itsuki Osawa, MD 
## Last updated on Jun 8th, 2023

## Load packages
library("tidyverse")
library("dplyr")
library("grf")

## Read the "Jseptic-DIC" cohort dataset
X <- Jseptic_DIC %>% dplyr::select(c("Age","Sex","APACHE2","SOFA_Lung","SOFA_Renal","SOFA_Coag","SOFA_liver","SOFA_CV",
                                     "SOFA_CNS","SOFA","WBC","Plt","PTINR","Fbg","FDP","Lac"))
Y <- as.vector(Jseptic_DIC$Survival_28day)
W <- as.vector(Jseptic_DIC$Treatment_PMX)

## Develop the machine learning-based causal forest model with 10-fold cross-fitting
## Please also see the following tutorial: "https://bookdown.org/halflearned/ml-ci-tutorial/hte-i-binary-treatment.html#data-driven-hypotheses"
num.rankings <- 5
num.folds <- 10
folds <- sort(seq(nrow(X)) %% num.folds) + 1

set.seed(1) 
CF <- causal_forest(X, Y, W, tune.parameters = "all", clusters = folds)

# Check the developed causal forest model
CF
test_calibration(CF)

# Compute AIPW scores based on the developed model
tau.hat <- predict(CF)$predictions # ITE (Figure 2A)
e.hat <- CF$W.hat # P[W=1|X]
m.hat <- CF$Y.hat # E[Y|X]

# Rank observations *within each fold into quintiles according to their CATE predictions.
ranking <- rep(NA, nrow(X))
for (fold in seq(num.folds)) {
  tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest=TRUE, labels=seq(num.rankings))
}

# Estimate mu.hat(X, 1) and mu.hat(X, 0) for obs in held-out sample
mu.hat.0 <- m.hat - e.hat * tau.hat        
mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat  

# Estimate AIPW scores and CATE in each ITE quintile group
library("lmtest")
library("sandwich")
aipw.scores <- tau.hat + W / e.hat * (Y -  mu.hat.1) - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0)
ols <- lm(aipw.scores ~ 0 + factor(ranking))
forest.ate <- data.frame("aipw", paste0("Quintile", seq(num.rankings)), coeftest(ols, vcov=vcovHC(ols, "HC2"))[,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL 
forest.ate ## Each CATE in each ITE quintile group (Figure 2B)

## Describe each quintile patient characteristics and identify determinants of ITEs (Table 1)
df_table <- cbind(Jseptic_DIC, ranking)
df_table %>% tbl_summary(by=ranking) %>% add_p() %>% bold_p() 

# Policytree (Supplemental Figure 3)
library("policytree")
dr.scores <- double_robust_scores(CF)
Potential_determinants <- Jseptic_DIC %>% dplyr::select(c("APACHE2","Plt","PTINR","Lac"))

Policy_tree <- policy_tree(Potential_determinants, dr.scores, depth = 2, split.step = 1)
Policy_tree