library(tidyverse)
library(caret)
library(matlib)

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

get_model_diagnostics <- function(model, moldesc) {
  
  # Model Diagnostics
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
  test_mols <- moldesc %>%
    filter(id %in% test_rows) %>%
    select(!c(id))
  pred_test <- predict(model, newdata = test_mols, se.fit = T)
  yhat_test <- pred_test$fit
  test_stdres <- pred_test$se.fit
  y_test = test_mols$pIC50
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
  
  y <- moldesc$pIC50
  y_pred <- predict(model, newdata = moldesc)
  
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
           mae_train = mae_train,
           mae_test = mae_test, 
           r_sq_test = r_sq_test,
           tropsha_3 = tropsha_3,
           tropsha_4 = tropsha_4,
           tropsha_5 = tropsha_5
  )
  )
}

desc_types <- c("blue","chemopy", "padel", "pybel", "rdkit")

SMLR_moldescs <- lapply(desc_types, FUN = function(desc_type) {
  desc_file <- paste0("SMLR/", desc_type, "_SMLR.csv")
  moldesc <- read_csv(desc_file)
})
names(SMLR_moldescs) <- desc_types

models <- lapply(desc_types, FUN = function(desc_type) {
  moldesc <- SMLR_moldescs[[desc_type]]
  
  # Fetch the SMLR-processed molecular descriptors with most sig variables
  train_mols <- moldesc %>%
    filter(id %in% train_rows) %>%
    select(!c(id))

  # Step-wise regression model
  lm <- lm(pIC50 ~ ., data = train_mols)
  slm <- step(lm)
  summary(slm)
  
  model <- slm
})
names(models) <- desc_types

model_diags <- lapply(desc_types, FUN = function(desc_type) {
  model <- models[[desc_type]]
  moldesc <- SMLR_moldescs[[desc_type]]
  get_model_diagnostics(model, moldesc)
})
names(model_diags) <- desc_types

# Hotelling's Test and Leverage (h)
# Pg.  37 (GUIDANCE DOCUMENT ON THE VALIDATION OF (QUANTITATIVE)STRUCTURE-ACTIVITY RELATIONSHIPS [(Q)SAR] MODELS)
Hotelling_AD <- function(alpha, train_moldesc, blinded_moldesc) {
  library(Hotelling)
  fit <- hotelling.test(train_moldesc, blinded_moldesc)
  if (fit$pval > alpha) {
    print("H0: Blinded and Train mol samples have the same mean")
  } else {
    print("H1: Blinded and Train mol samples have different means")
  }
  return (fit$pval > alpha)
}

# Leverage and Williams Plot
Leverage_AD <- function(model, train_moldesc, test_moldesc, blinded_moldesc) {
  X <- as.matrix(train_moldesc)
  X_sq_inv <- Inverse(t(X) %*% X)
  
  # Warning Leverage
  h_warn <- 3*(ncol(X) + 1)/nrow(X) 
  
  # Williams Plot
  desc <- train_moldesc
  h_train <- apply(desc, 1, FUN = function(x) {
    t(x) %*% X_sq_inv %*% x
  })
  
  train_stdres <- rstandard(model)
  w_train <- data.frame(stdres = train_stdres, h = h_train, class = "train")
  rownames(w_train) <- paste("T",seq_along(h_train), sep="")
  desc <- test_moldesc
  h_test <- apply(desc, 1, FUN = function(x) {
    t(x) %*% X_sq_inv %*% x
  })
  
  pred_test <- predict(model, newdata = test_moldesc, se.fit = T)
  test_stdres <- pred_test$se.fit
  w_test <- data.frame(stdres = test_stdres, h = h_test, class = "test")
  rownames(w_test) <- paste("S",seq_along(h_test), sep="")
  w <- rbind(w_train, w_test)
  
  print(ggplot(data = w, aes(x = h, y = stdres, shape = class, label = rownames(w))) +
    geom_point() +
    geom_text() +
    geom_hline(yintercept = mean(w_train$stdres) + 2*sd(w_train$stdres), color = "red") +
    geom_hline(yintercept = mean(w_train$stdres) - 2*sd(w_train$stdres), color = "red") +
    geom_vline(xintercept = h_warn, color = "red")
  )
  desc <- blinded_moldesc
  h <- apply(desc, 1, FUN = function(x) {
    t(x) %*% X_sq_inv %*% x
  })
  
  in_AD_Leverage_test = which(h_test <= h_warn)
  out_AD_Leverage_test = which(h_test > h_warn)
  
  return(list(in_AD_Leverage_test, out_AD_Leverage_test))
}

# Multivariate Chebyshev Inequality with Estimated Mean and Variance by Stellato
multivariate_chebychev <- function(train_test_mols, ext_mols) {
  n = ncol(train_test_mols)
  N = nrow(train_test_mols)
  
  mu = colMeans(train_test_mols)
  Sigma = cov(train_test_mols)
  
  ds_sq <- apply(ext_mols, 1, function(x) {mahalanobis(x, mu, Sigma)})
  probs <- sapply(ds_sq, FUN = function(lambda) {min(1, n*(N^2-1+N*lambda^2)/(N^2*lambda^2))})
  return(probs)
}

Chebychev_AD <- function(train_moldesc, test_moldesc, blinded_moldesc) {
  probs_train <- multivariate_chebychev(train_moldesc, train_moldesc)
  probs_test <- multivariate_chebychev(train_moldesc, test_moldesc)
  probs <- multivariate_chebychev(train_moldesc, blinded_moldesc)
  
  in_AD_Chebychev_test <- which(probs_test >= min(probs_train) & probs_test <= max(probs_train))
  out_AD_Chebychev_test <- seq_along(probs_test)[-in_AD_Chebychev_test]
  
  in_AD_Chebychev <- which(probs >= min(probs_train) & probs <= max(probs_train))
  out_AD_Chebychev <- seq_along(probs)[-in_AD_Chebychev]
  
  return(list(in_AD_Chebychev_test, out_AD_Chebychev_test, in_AD_Chebychev, out_AD_Chebychev))
}

desc_files <- c("blue_desc.csv", "chemopy_desc.csv", "padel_desc.csv", "pybel_desc.csv", "rdkit_desc.csv")
ADs <- lapply(desc_types, FUN = function(desc_type) {
  moldesc <- read_csv(paste0("descriptors/", desc_type, "_desc.csv")) 
  model <- models[[desc_type]]
  model_vars <- gsub('\`', '', variable.names(model)[-1])
  
  train_moldesc <- moldesc[train_rows, model_vars]
  test_moldesc <- moldesc[test_rows, model_vars]
  train_test_moldesc <- moldesc[train_test_rows, model_vars]
  blinded_moldesc <- moldesc[blinded_rows, model_vars]
  
  alpha = 0.05
  hotelling_AD <- Hotelling_AD(alpha, train_moldesc, blinded_moldesc)
  
  leverage_AD <- Leverage_AD(model, train_moldesc, test_moldesc, blinded_moldesc)
  
  chebychev_AD <- Chebychev_AD(train_moldesc, test_moldesc, blinded_moldesc)
  
  
  
})





