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


desc_types <- c("blue","chemopy", "padel", "pybel", "rdkit", "padel_rdkit", "all")

padel <- read_csv(paste0("descriptors/", "padel_desc.csv"))
rdkit <- read_csv(paste0("descriptors/", "rdkit_desc.csv"))
padel_rdkit <- cbind(padel, rdkit)
write_csv(padel_rdkit, "descriptors/padel_rdkit_desc.csv")

files <- paste0("descriptors/", desc_types[1:5], "_desc.csv")
all <- do.call("cbind", lapply(files, read_csv))
write_csv(all, "descriptors/all_desc.csv")

SMLR_path = "~/ddh/S-MLR 1.2_24March2017/"
descs <- lapply(desc_types, FUN = function(desc_type) {
  desc_file <- paste0(desc_type, "_desc.csv")
  df <- 
    read_csv(paste0("descriptors/", desc_file)) %>%
    select(where(is.numeric)) %>%
    mutate(id = 1:nrow(.)) %>%
    relocate(id) %>%
    mutate(pIC50 = smiles$pIC50[id]) %>%
    relocate(pIC50, .after = last_col()) %>%
    filter(!is.na(pIC50))
  write_csv(df, path = paste0(SMLR_path, "data/", desc_file))
  df
})
names(descs) <- desc_types

# TODO: Run SMLR tool programmatically 
SMLR_moldescs <- lapply(desc_types, FUN = function(desc_type) {
  desc_file <- paste0(SMLR_path, "output/", desc_type, "_SMLR.csv")
  moldesc <- read_csv(desc_file)
})
names(SMLR_moldescs) <- desc_types

scalers <- lapply(desc_types, FUN = function(desc_type) {
  SMLR_moldescs[[desc_type]] %>%
    select(!c(id, pIC50)) %>%
    preProcess(method = "range", rangeBounds = c(-1,1))
})
names(scalers) <- desc_types

SMLR_models <- lapply(desc_types, FUN = function(desc_type) {
  moldesc <- SMLR_moldescs[[desc_type]]
  scaler <- scalers[[desc_type]]
  
  # Fetch the SMLR-processed molecular descriptors with most sig variables
  train_mols <- moldesc %>%
    filter(id %in% train_rows) %>%
    select(!c(id)) %>%
    predict(scaler, newdata = .)

  # Step-wise regression model
  lm <- lm(pIC50 ~ ., data = train_mols, na.action = "na.fail")
  slm <- step(lm)
  summary(slm)
  
  model <- slm
})
names(SMLR_models) <- desc_types

get_model_diagnostics <- function(desc_type) {
  model <- SMLR_models[[desc_type]]
  moldesc <- SMLR_moldescs[[desc_type]]
  scaler <- scalers[[desc_type]]
  
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
    select(!c(id)) %>%
    predict(scaler, newdata = .)
  
  pred_test <- predict(model, newdata = test_mols, se.fit = T)
  yhat_test <- pred_test$fit
  test_stdres <- pred_test$se.fit
  y_test = test_mols$pIC50
  cor(yhat_test, y_test)
  plot(y_test, yhat_test)
  test_obs_pred <- data.frame(y_test = y_test, yhat_test = yhat_test)
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
  y_pred <- predict(model, newdata = predict(scaler, moldesc))
  
  tropsha_k = sum(y * y_pred)/sum(y_pred^2)
  tropsha_k_dash = sum(y * y_pred)/sum(y^2)
  
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
  print(paste("3b: k =", tropsha_k, ", 0.85 <= k <= 1.15:", tropsha_k >= 0.85 & tropsha_k <= 1.15))
  print(paste("4a: (r\'_sq-r0\'_sq)/r\'_sq =", tropsha_4, "< 0.1:", tropsha_4 < 0.1))
  print(paste("4b: k\' =", tropsha_k_dash, ", 0.85 <= k\' <= 1.15:", tropsha_k_dash >= 0.85 & tropsha_k_dash <= 1.15))
  print(paste("5: |r0_sq - r0\'_sq| =", tropsha_5, " < 0.3:", tropsha_5 < 0.3))
  
  tropsha_cond <- tropsha_3 < 0.1 &&
    (tropsha_k >= 0.85 & tropsha_k <= 1.15) &&
    tropsha_4 < 0.1 &
    (tropsha_k_dash >= 0.85 && tropsha_k_dash <= 1.15) &&
    tropsha_5 < 0.3
  
  fit_conds <- rsq_train > 0.7 &&
    q_sq_train > 0.7 &&
    q_sq_ext > 0.7 &&
    r_sq_test > 0.6
  
  return(c(rsq_train = rsq_train,
           see_train = see_train,
           q_sq_train = q_sq_train,
           test_obs_pred = test_obs_pred,
           q_sq_ext = q_sq_ext,
           mae_train = mae_train,
           mae_test = mae_test, 
           r_sq_test = r_sq_test,
           tropsha_3 = tropsha_3,
           tropsha_k = tropsha_k,
           tropsha_4 = tropsha_4,
           tropsha_k_dash = tropsha_k_dash,
           tropsha_5 = tropsha_5,
           tropsha_cond = tropsha_cond,
           fit_conds = fit_conds
  )
  )
}

model_diags <- lapply(desc_types, FUN = function(desc_type) {
  get_model_diagnostics(desc_type)
})
names(model_diags) <- desc_types

# XternalValidationPlus tool - for MAE-based criteria tool
library(readxl)
XVP_path <- "~/ddh/XternalValidationPlus_Updated13April16/"
y_train <- smiles$pIC50[train_rows]
mu_y_train <- mean(y_train)
mu_y_train
range_y_train <- diff(range(y_train))
range_y_train

mae_crits <- sapply(desc_types, FUN = function(desc_type) {
  print(desc_type)
  diag <- model_diags[[desc_type]]
  df <- data.frame(Yobs_test = diag$test_obs_pred.y_test,
                   Ypred_test = diag$test_obs_pred.yhat_test)
  write_excel_csv(df, paste0(XVP_path, "/data/ytest_", desc_type, ".csv"))
  
  # Run the jar manually
  system(paste("java -Xmx4000m -jar", paste0(XVP_path, "XternalValidationPlusUpdated.jar")))
  
  df <- read_excel(paste0(XVP_path, "/output/", desc_type, "_ExternalValidation.xls"))
  mae_crit <- toupper(tail(df,1)[3])
})


# Select Desc Types based on DDH01 problem statement:
# R2 > 0.7, LOO-Q2 > 0.7, Q2ext_F1 > 0.7, Tropsha’s criteria: Pass, 
# MAE based criteria: Moderate or good
diag_flags <- sapply(desc_types, FUN = function(desc_type) {
  diag <- model_diags[[desc_type]]
  
  mae_crit <- mae_crits[desc_type]
  mae_cond <- mae_crit != 'BAD'
  
  diag[["fit_conds"]] && diag[["tropsha_cond"]] && mae_cond
    
})

# Choose the desc types that passed all the criteria
sel_desc_types <- desc_types[diag_flags]

# Run a random forest on the predictions from the selected models
library(caret)
get_ensemble_model <- function(sel_desc_types, moldescs, models) {
  y_hats <- cbind(sapply(sel_desc_types, FUN = function(desc_type) {
    moldesc <- moldescs[[desc_type]]
    model <- models[[desc_type]]
    scaler <- scalers[[desc_type]]
    
    predict(model, newdata = predict(scaler, moldesc))
  }))
  y_hats <- data.frame(y_hats)
  print(y_hats)
  y_hats$pIC50 <- smiles$pIC50[train_test_rows]
  
  rf_fit <- train(pIC50 ~ ., data = y_hats, method = "rf", metric = "RMSE")
}

# Train the ensemble
rf_fit <- get_ensemble_model(sel_desc_types, SMLR_moldescs, SMLR_models)
y_pred <- predict(rf_fit)



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
Leverage_AD <- function(model, train_moldesc, test_moldesc, type, plot_train = TRUE) {
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
  w_test <- data.frame(stdres = test_stdres, h = h_test, class = type)
  test_obs_type = switch(type,
                         test = 'S',
                         blinded = 'B',
                         cas = 'C')
  rownames(w_test) <- paste(test_obs_type,seq_along(h_test), sep="")
  if (plot_train) {
    w <- rbind(w_train, w_test)
  } else {
    w <- w_test
  }
  
  print(ggplot(data = w, aes(x = h, y = stdres, shape = class, label = rownames(w))) +
    geom_point() +
    geom_text() +
    geom_hline(yintercept = mean(w_train$stdres) + 2*sd(w_train$stdres), color = "red") +
    geom_hline(yintercept = mean(w_train$stdres) - 2*sd(w_train$stdres), color = "red") +
    geom_vline(xintercept = h_warn, color = "red")
  )

  in_AD_Leverage = which(h_test <= h_warn)
  out_AD_Leverage = which(h_test > h_warn)
  
  return(list(in_AD_Leverage, out_AD_Leverage))
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

Chebychev_AD <- function(train_moldesc, test_moldesc) {
  probs_train <- multivariate_chebychev(train_moldesc, train_moldesc)
  probs_test <- multivariate_chebychev(train_moldesc, test_moldesc)
  
  in_AD_Chebychev <- which(probs_test >= min(probs_train) & probs_test <= max(probs_train))
  out_AD_Chebychev <- seq_along(probs_test)[-in_AD_Chebychev]
  
  return(list(in_AD_Chebychev, out_AD_Chebychev))
}

Tanimoto_AD <- function() {
  # Tanimoto AD
  
  train_smiles <- parse.smiles(smiles$SMILES[train_rows])
  test_smiles <- parse.smiles(smiles$SMILES[test_rows])
  train_test_smiles <- parse.smiles(smiles$SMILES[train_test_rows])
  blinded_smiles <- parse.smiles(smiles$SMILES[blinded_rows])
  
  # 1-NN in Training based on Test pIC50 (activity)
  y_train <- smiles$pIC50[train_rows]
  y_test <- smiles$pIC50[test_rows]
  activ_1nn <- sapply(y_test, FUN = function(y) {which.min(abs(y_train - y))})
  activ_1nn_sim <- sapply(seq_along(activ_1nn), FUN = function(i) {
    y <- y_test[i]
    row <- activ_1nn[i]
    abs(y_train[row] - y)
  })
  
  # Best sig_type using RMSE between 1-NN Activity and 1-NN Tanimoto similarity on training set
  circular_types = c('ECFP0','ECFP2','ECFP4','ECFP6','FCFP0','FCFP2','FCFP4','FCFP6')
  sig_types = c('standard', 'extended', 'graph', 'hybridization', 
                'maccs', 'estate', 'pubchem', 'shortestpath', 
                'substructure', paste('circular', circular_types, sep = '.'))
  options("java.parameters"=c("-Xmx4000m"))
  rmse <- sapply(sig_types, FUN = function(sig_type) {
    
    if (str_starts(sig_type, 'circular')) {
      circ_type = strsplit(sig_type, '.', fixed = T)[[1]][2]
      fps_train <- lapply(train_smiles, get.fingerprint, type='circular', circular.type=circ_type)
      fps_test <- lapply(test_smiles, get.fingerprint, type='circular', circular.type=circ_type)
    } else {
      fps_train <- lapply(train_smiles, get.fingerprint, type=sig_type)
      fps_test <- lapply(test_smiles, get.fingerprint, type=sig_type)
    }
    
    test_sim <- fingerprint::fp.sim.matrix(fps_train, fps_test, method='tanimoto')
    
    # 1-NN Tanimoto
    tanimoto_1nn <- apply(test_sim, 2, which.max)
    tanimoto_1nn_sim <- sapply(seq_along(tanimoto_1nn_test), FUN = function(col) {
      test_sim[tanimoto_1nn[col], col]
    })
    
    # Find the activity similarity based on 1-NN Tanimoto
    tanimoto_1nn_activ_sim <- sapply(seq_along(tanimoto_1nn), FUN = function(i) {
      y <- y_test[i]
      row <- tanimoto_1nn[i]
      abs(y_train[row] - y)
    })
    
    # RMSE
    # sqrt(mean(tanimoto_1nn_activ_sim^2))
    cor(tanimoto_1nn_activ_sim, tanimoto_1nn_sim)
    #sqrt(mean((activ_1nn_sim - tanimoto_1nn_activ_sim)^2))
  })
  
  #sig_type = names(which.min(rmse))
  sig_type = 'maccs'
  if (str_starts(sig_type, 'circular')) {
    circ_type = strsplit(sig_type, '.', fixed = T)[[1]][2]
    fps_train <- lapply(train_smiles, get.fingerprint, type='circular', circular.type=circ_type)
    fps_test <- lapply(test_smiles, get.fingerprint, type='circular', circular.type=circ_type)
    fps_blinded <- lapply(blinded_smiles, get.fingerprint, type='circular', circular.type=circ_type)
  } else {
    fps_train <- lapply(train_smiles, get.fingerprint, type=sig_type)
    fps_test <- lapply(test_smiles, get.fingerprint, type=sig_type)
    fps_blinded <- lapply(blinded_smiles, get.fingerprint, type=sig_type)
  }
  
  test_sim <- fingerprint::fp.sim.matrix(fps_train, fps_test, method='tanimoto')
  tanimoto_1nn_test <- apply(test_sim, 2, which.max)
  tanimoto_1nn_sim_test <- sapply(seq_along(tanimoto_1nn_test), FUN = function(col) {
    test_sim[tanimoto_1nn_test[col], col]
  })
  
  blinded_sim <- fingerprint::fp.sim.matrix(fps_train, fps_blinded, method='tanimoto')
  tanimoto_1nn <- apply(blinded_sim, 2, which.max)
  tanimoto_1nn_sim <- sapply(seq_along(tanimoto_1nn), FUN = function(col) {
    blinded_sim[tanimoto_1nn[col], col]
  })
  
  sim_threshold <- 0.8
  in_AD_Tanimoto_test <- which(tanimoto_1nn_sim_test >= sim_threshold)
  out_AD_Tanimoto_test <- which(tanimoto_1nn_sim_test < sim_threshold)
  
  in_AD_Tanimoto <- which(tanimoto_1nn_sim >= sim_threshold)
  out_AD_Tanimoto <- which(tanimoto_1nn_sim < sim_threshold)
  
  return(list(in_AD_Tanimoto_test, out_AD_Tanimoto_test, in_AD_Tanimoto, out_AD_Tanimoto))
}

ADS_path <- "~/ddh/ADUsingStdApproach/"
Standardization_AD <- function(desc_type, train_moldesc, test_moldesc, type = "test") {
  print(desc_type)
  train_moldesc %>%
    mutate(id = 1: nrow(.)) %>%
    relocate(id) %>%
    write_excel_csv(path = paste0(ADS_path, "data/", desc_type, "_train.csv"), na = "")
  
  test_moldesc %>%
    mutate(id = 1: nrow(.)) %>%
    relocate(id) %>%
    write_excel_csv(path = paste0(ADS_path, "data/", desc_type, "_", type, ".csv"), na ="")
  
  # Run the jar manually
  system(paste("java -Xmx4000m -jar", paste0(ADS_path, "ADUsingStdApproach.jar")))
  
  AD_Standardization <- 
    read_csv(paste0(ADS_path, "output/", desc_type,"_", type, "_Test_AD.csv")) %>%
    select(last_col()) 
  
  out_AD_Standardization = which(AD_Standardization == "Outside AD")
  in_AD_Standardization = which(AD_Standardization != "Outside AD")
  
  return(list(in_AD_Standardization, out_AD_Standardization))
  
}

# Finally !
sel_desc_types <- desc_types[diag_flags]

AD_preds <- lapply(sel_desc_types, FUN = function(desc_type) {
  moldesc <- read_csv(paste0("descriptors/", desc_type, "_desc.csv")) 
  model <- SMLR_models[[desc_type]]
  model_vars <- gsub('\`', '', variable.names(model)[-1])
  scaler <- scalers[[desc_type]]
  
  train_moldesc <- 
    moldesc[train_rows,] %>% 
    predict(scaler, newdata = .) %>%
    select(all_of(model_vars))
  test_moldesc <- 
    moldesc[test_rows,] %>% 
    predict(scaler, newdata = .) %>%
    select(all_of(model_vars))
  blinded_moldesc <- 
    moldesc[blinded_rows,] %>% 
    predict(scaler, newdata = .) %>%
    select(all_of(model_vars))
  
  alpha = 0.05
  hotelling_AD <- Hotelling_AD(alpha, train_moldesc, blinded_moldesc)
  
  AD_Leverage_test <- Leverage_AD(model, train_moldesc, test_moldesc, type = "test")
  AD_Leverage_blinded <- Leverage_AD(model, train_moldesc, blinded_moldesc, type = "blinded")
  
  AD_Chebychev_test <- Chebychev_AD(train_moldesc, test_moldesc)
  AD_Chebychev_blinded <- Chebychev_AD(train_moldesc, blinded_moldesc)
  
  AD_Standardization_test <- Standardization_AD(desc_type, 
                                           train_moldesc, 
                                           test_moldesc,
                                           type = "test")
  AD_Standardization_blinded <- Standardization_AD(desc_type, 
                                                train_moldesc, 
                                                blinded_moldesc,
                                                type = "blinded")
  
  pred_blinded <- predict(model, newdata = blinded_moldesc)

  return(list(pred_blinded = pred_blinded,
              AD_Leverage_test = AD_Leverage_test,
              AD_Leverage_blinded = AD_Leverage_blinded,
              AD_Chebychev_test = AD_Chebychev_test,
              AD_Chebychev_blinded = AD_Chebychev_blinded,
              AD_Standardization_test = AD_Standardization_test,
              AD_Standardization_blinded = AD_Standardization_blinded
              )
         )
  
})
names(AD_preds) <- sel_desc_types

# Ensemble Blinded Predictions and AD Info
ensemble_predict <- function(models, desc_files) {
  preds <- mapply(FUN = function(model, desc_file) {
    read_csv(desc_file) %>%
      predict(model, newdata = .)
  }, models, desc_files)
  predict(rf_fit, newdata = preds)
}

#preds_blinded <- mapply(predict, models[sel_desc_types], descs[sel_desc_types][blinded_rows])
preds_blinded <- cbind(sapply(sel_desc_types, FUN = function(desc_type) {
  AD_preds[[desc_type]][["pred_blinded"]]
}))
pred_blinded_ensemble <- predict(rf_fit, newdata = preds_blinded) 

# Ensemble Blinded AD Information
out_AD_blinded <- Reduce("intersect", lapply(sel_desc_types, FUN = function(desc_type) {
  AD <- AD_preds[[desc_type]]
  out_AD_Leverage_blinded <- AD[["AD_Leverage_blinded"]][[2]]
  out_AD_Chebychev_blinded <- AD[["AD_Chebychev_blinded"]][[2]]
  out_AD_Standardization_blinded <- AD[["AD_Standardization_blinded"]][[2]]
  
  Reduce("union", list(out_AD_Leverage_blinded,
                       out_AD_Chebychev_blinded,
                       out_AD_Standardization_blinded
  ))
}))

# Ensemble Test AD Information
out_AD_test <- Reduce("intersect", lapply(sel_desc_types, FUN = function(desc_type) {
  AD <- AD_preds[[desc_type]]
  out_AD_Leverage_test <- AD[["AD_Leverage_test"]][[2]]
  out_AD_Chebychev_test <- AD[["AD_Chebychev_test"]][[2]]
  out_AD_Standardization_test <- AD[["AD_Standardization_test"]][[2]]
  
  Reduce("union", list(out_AD_Leverage_test,
                       out_AD_Chebychev_test,
                       out_AD_Standardization_test
  ))
}))

# Write predictions and AD information to CSV files

# DDT Input Form 2
preds_df <- 
  data.frame(id = blinded_rows, pIC50_pred = pred_blinded_ensemble) %>%
  mutate(AD = if_else(id %in% blinded_rows[out_AD_blinded], "Outside AD", "-")) %>%
  write_excel_csv("results/Input Form 2_DDT1-01.csv")

# DDT Input Form 1
sapply(sel_desc_types, FUN = function(desc_type) {
  model <- SMLR_models[[desc_type]]
  model_vars <- names(coef(model))[-1]
  diag <- model_diags[[desc_type]]
  
  csv_file <- paste0("results/Input Form 1_DDT1-01_", desc_type, ".csv")
  
  df <- 
    data.frame(diag = c("DDT1-01",
                        diag$rsq_train,
                        diag$see_train,
                        diag$q_sq_train,
                        diag$mae_train,
                        diag$q_sq_ext,
                        diag$mae_test,
                        ifelse(diag$tropsha_cond, "Yes", "No"),
                        mae_crits[desc_type],
                        paste(train_rows[out_AD_test])
                        )
    ) %>%
    write_excel_csv(csv_file)
  
  # Insert empty row
  write_excel_csv(data.frame(c(NA, "Training Set")), csv_file, na = "", append = T)
  
  moldesc <- SMLR_moldescs[[desc_type]] 
  train_moldesc <- 
    moldesc %>% 
    filter(id %in% train_rows) %>%
    select(append(model_vars, c("id", "pIC50"))) %>%
    relocate(id) %>%
    relocate(pIC50, .after = last_col()) %>%
    mutate(pIC50pred = predict(model, newdata = .)) %>%
    write_excel_csv(csv_file, append = T, col_names = T)
  
  # Insert empty row
  write_excel_csv(data.frame(c(NA, "Test Set")), csv_file, na = "", append = T)
  
  test_moldesc <- 
    moldesc %>% 
    filter(id %in% test_rows) %>%
    select(append(model_vars, c("id", "pIC50"))) %>%
    relocate(id) %>%
    relocate(pIC50, .after = last_col()) %>%
    mutate(pIC50pred = predict(model, newdata = .)) %>%
    write_excel_csv(csv_file, append = T, col_names = T)
  
  
})


# Sanity check using MACCS molecular fingerprint
library(rcdk)
train_smiles <- parse.smiles(smiles$SMILES[train_rows])
test_smiles <- parse.smiles(smiles$SMILES[test_rows])
train_test_smiles <- parse.smiles(smiles$SMILES[train_test_rows])
blinded_smiles <- parse.smiles(smiles$SMILES[blinded_rows])

sig_type = 'maccs'
fps_train <- lapply(train_smiles, get.fingerprint, type=sig_type)
fps_test <- lapply(test_smiles, get.fingerprint, type=sig_type)
fps_blinded <- lapply(blinded_smiles, get.fingerprint, type=sig_type)

blinded_sim <- fingerprint::fp.sim.matrix(fps_train, fps_blinded, method='tanimoto')
tanimoto_1nn <- apply(blinded_sim, 2, which.max)
pIC50_1nn <- smiles$pIC50[train_rows[tanimoto_1nn]]
tanimoto_1nn_sim <- sapply(seq_along(tanimoto_1nn), FUN = function(col) {
  blinded_sim[tanimoto_1nn[col], col]
})
df <- data.frame(yhat = pred_blinded_ensemble, y = pIC50_1nn, sim = tanimoto_1nn_sim)
df$yerr <- df$y - df$yhat

ggplot(df, aes(yhat, y)) +
  geom_point(aes(size = sim)) +
  geom_label(aes(label = paste0('B', 1:nrow(df))), nudge_y = 0.25)




