library(tidyverse)
library(caret)

smiles <- read_csv("smiles.csv") %>%
  mutate(pIC50 = ifelse(pIC50 == "BLINDED", NA, pIC50)) %>%
  mutate(pIC50 = as.numeric(pIC50))

all_rows = 1:nrow(smiles)
train_test_rows = which(!is.na(smiles$pIC50))#[-c(32,43)]
blinded_rows = all_rows[is.na(smiles$pIC50)]

# Train-Test split
set.seed(1024)
test_rows = sample(train_test_rows, 0.3*floor(length(train_test_rows)))
train_rows = setdiff(train_test_rows, test_rows)


desc_types <- c("blue_desc","chemopy", "padel", "pybel", "rdkit")

models <- lapply(desc_types, FUN = function(desc_type) {
  desc_file <- paste0("SMLR/", desc_type, "_SMLR.csv")
  do_regression(desc_file)
})
names(models) <- desc_types

do_regression <- function(desc_file) {
  # Fetch the SMLR-processed molecular descriptors with most sig variables
  mols <- read_csv(desc_file)
  train_mols <- mols %>%
                  filter(id %in% train_rows) %>%
                    select(!c(id))
                    
  
  # Step-wise regression model
  lm <- lm(pIC50 ~ ., data = train_mols)
  slm <- step(lm)
  summary(slm)
  
  model <- slm
  return(model)
}

get_model_diagnostics <- function(model) {
  # Model Diagnostics
  model <- slm
  pval = 1
  coefs <- coef(summary(model))
  model_vars <- names(which(coefs[-1,4] <= pval))
  plot(model)
  
  # R-squared (Training)
  rsq_train = summary(model)$r.squared
  print(paste("R-squared (Training) = ", rsq_train))
  
  # Standard Error of Estimation (Training)
  y_train <- smiles$pIC50[train_rows]
  pred_train <- predict(model, se.fit = T)
  yhat_train <- pred_train$fit
  see_train <- sqrt(mean((y_train - yhat_train)^2))
  print(paste("Standard Error of Estimation (Training) = ", see_train))
  
  # Empirical CDF/Density and Shapiro-Wilkes test for normality of residuals
  train_stdres = rstandard(model)
  qqnorm(train_stdres, 
         ylab="Standardized Residuals", 
         xlab="Normal Scores", 
         main="DDH-01-Padel") 
  qqline(train_stdres)
  
  # Shapiro test for normality of residuals
  plot(ecdf(train_stdres))
  plot(density(train_stdres, na.rm = T))
  shapiro.test(train_stdres)
  
  # LOO-CV residuals from rstandard help page
  err_cv = rstandard(model, type="predictive")^2
  PRESS_train <- sum(err_cv, na.rm = T)
  tot_sq = sum((y_train - mean(y_train))^2)
  q_sq_train = 1 - PRESS_train/tot_sq
  print(paste("LOO-Q2 (Training) = ", q_sq_train))
  
  # Q2 ext F1
  pred_test <- predict(model, newdata = moldesc[test_rows, model_vars], se.fit = T)
  yhat_test <- pred_test$fit
  test_stdres <- pred_test$se.fit
  y_test = smiles$pIC50[test_rows]
  cor(yhat_test, y_test)
  plot(y_test, yhat_test)
  pred_test <- data.frame(y_test = y_test, yhat_test = yhat_test)
  q_sq_ext = 1.0 - sum((y_test - yhat_test)^2)/sum((y_test - mean(y_train))^2)
  print(paste("Q2extF1 = ", q_sq_ext))
  
  # MAE and MAPE
  mae_train <- mean(abs(y_train - yhat_train))
  print(paste("MAE (Training) = ", mae_train))
  mae_test <- mean(abs(y_test - yhat_test))
  print(paste("MAE (Test) = ", mae_test))
  mape <- mean(abs((y_test - yhat_test)/y_test))
  mape
  
  # Tropsha's Criteria
  r_sq_test = 1 - sum((y_test - yhat_test)^2)/sum((y_test - mean(y_test))^2)
  
  y <- smiles$pIC50[train_test_rows]
  y_pred <- predict(model, newdata = moldesc[train_test_rows, model_vars])
  
  k = sum(y * y_pred)/sum(y_pred^2)
  k_dash = sum(y * y_pred)/sum(y^2)
  
  lm_obs_pred <- lm(y ~ y_pred)
  summary(lm_obs_pred)
  r_sq = summary(lm_obs_pred)$r.squared
  lm_obs_pred0 <- lm(y ~ y_pred - 1)
  summary(lm_obs_pred0)
  r0_sq = summary(lm_obs_pred0)$r.squared
  tropsha_3 = (r_sq-r0_sq)/r_sq
  
  lm_pred_obs <- lm(y_pred ~ y)
  summary(lm_pred_obs)
  r_dash_sq = summary(lm_pred_obs)$r.squared
  lm_pred_obs0 <- lm(y_pred ~ y - 1)
  summary(lm_pred_obs0)
  r0_dash_sq = summary(lm_pred_obs0)$r.squared
  tropsha_4 = (r_dash_sq-r0_dash_sq)/r_dash_sq
  
  tropsha_5 = abs(r0_sq - r0_dash_sq)
  
  print("Tropsha's Conditions:")
  print(paste("1: Q-squared (Training)=", q_sq_train,"> 0.5", q_sq_train > 0.5))
  print(paste("2: R-squared (Test)=", r_sq_test,"> 0.6", r_sq_test > 0.6))
  print(paste("3a: (r_sq-r0_sq)/r_sq =", tropsha_3, "< 0.1:", tropsha_3 < 0.1))
  print(paste("3b: k =", k, ", 0.85 <= k <= 1.15:", k >= 0.85 & k <= 1.15))
  print(paste("4a: (r\'_sq-r0\'_sq)/r\'_sq =", tropsha_4, "< 0.1:", tropsha_4 < 0.1))
  print(paste("4b: k\' =", k_dash, ", 0.85 <= k\' <= 1.15:", k_dash >= 0.85 & k_dash <= 1.15))
  print(paste("5: |r0_sq - r0\'_sq| =", tropsha_5, " < 0.3:", tropsha_5 < 0.3))
  
  
  return(c(rsq_train = rsq_train,
           see_train = see_train,
           q_sq_train = q_sq_train,
           pred_test = pred_test,
           q_sq_ext = q_sq_ext,
           mape = mape,
           r_sq_test = r_sq_test,
           tropsha_3 = tropsha_3,
           tropsha_4 = tropsha_4,
           tropsha_5 = tropsha_5
           )
         )
}









