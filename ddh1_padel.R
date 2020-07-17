library(tidyverse)
library(caret)

smiles <- read_csv("smiles.csv") %>%
  mutate(pIC50 = ifelse(pIC50 == "BLINDED", NA, pIC50)) %>%
  mutate(pIC50 = as.numeric(pIC50))

padel <- read_csv("padel_descriptors_1d_2d.csv")
moldesc <- padel[, -1]
vars <- names(moldesc)

mols <- cbind(smiles, moldesc)

all_rows = 1:nrow(smiles)
train_test_rows = which(!is.na(smiles$pIC50))#[-c(32,43)]
blinded_rows = all_rows[is.na(smiles$pIC50)]

# Activity correlation
activ_cors = apply(moldesc[train_test_rows,], 2, FUN = function(v) {
  r = 0
  tryCatch({
    r = cor(smiles$pIC50[train_test_rows], v, use = "pairwise.complete.obs")
  }, finally = {return(r)})
}) 

useSMLR = T
if (useSMLR) {
  df <- read_csv("output_SMLR.csv")
  sel_vars <- names(df)[names(df) %in% vars]
} else {
  # Select top-n vars based on correlation with activity from padel dataset
  nvars = ceiling(length(train_test_rows)/1.25) * 0.7
  sel_vars = vars[head(sort.list(abs(activ_cors), decreasing = T), nvars)]
  
  # Cross Correlation
  x_cors = cor(moldesc[train_test_rows, sel_vars])
  hi_cors <- which(x_cors > 0.5, arr.ind = T)
}

  
# Train-Test split
set.seed(1024)
test_rows = sample(train_test_rows, 0.3*floor(length(train_test_rows)))
train_rows = setdiff(train_test_rows, test_rows)

# Prepare data for Regression
train_mols = mols[train_rows, append(sel_vars, "pIC50")]

# Step-wise regression model
lm <- lm(pIC50 ~ ., data = train_mols)
slm <- step(lm)
summary(slm)

# Model Diagnostics
model <- slm
pval = 1
coefs <- coef(summary(model))
model_vars <- names(which(coefs[-1,4] <= pval))
plot(model)

# R-squared and adjusted R-squared
rsq = summary(model)$r.squared
print(paste("R2 = ", rsq))

# Empirical CDF/Density and Shapiro-Wilkes test for normality of residuals
train_stdres = rstandard(model)
qqnorm(train_stdres, 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="DDH-01-Padel") 
qqline(train_stdres)

plot(ecdf(stdres))
plot(density(stdres, na.rm = T))
shapiro.test(stdres)

# LOO-CV residuals from rstandard help page
err_cv = rstandard(model, type="predictive")^2
PRESS <- sum(err_cv, na.rm = T)
y_cv = y[!is.nan(err_cv)]
tot_sq = sum((y_cv - mean(y_train))^2)
q_sq = 1 - PRESS/tot_sq
print(paste("LOO-Q2 = ", q_sq))

# Q2 ext F1
pred_test <- predict(model, newdata = moldesc[test_rows, model_vars], se.fit = T)
yhat_test <- pred_test$fit
test_stdres <- pred_test$se.fit
y_test = smiles$pIC50[test_rows]
cor(yhat_test, y_test)
plot(y_test, yhat_test)
pred_test <- data.frame(y_test = y_test, yhat_test = yhat_test)
write_csv(pred_test, "pred_test.csv")
y_train = smiles$pIC50[train_rows]
q_sq_ext = 1.0 - sum((y_test - yhat_test)^2)/sum((y_test - mean(y_train))^2)
print(paste("Q2extF1 = ", q_sq_ext))

# MAE and MAPE
mae <- mean(abs(y_test - yhat_test))
print(paste("MAE = ", mae))
mape <- mean(abs((y_test - yhat_test)/y_test))
mape

# Tropsha's Criteria
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
print(paste("Tropsha Condition 3 = ", (r_sq-r0_sq)/r_sq, " and k = ", k))

lm_pred_obs <- lm(y_pred ~ y)
summary(lm_pred_obs)
r_dash_sq = summary(lm_pred_obs)$r.squared
lm_pred_obs0 <- lm(y_pred ~ y - 1)
summary(lm_pred_obs0)
r0_dash_sq = summary(lm_pred_obs0)$r.squared
print(paste("Tropsha Condition 4 = ", (r_dash_sq-r0_dash_sq)/r_dash_sq, " and k\' = ", k_dash))

print(paste("Tropsha Condition 5 = ", abs(r0_sq - r0_dash_sq)))

# Final Variables
pval = 1
coefs <- coef(summary(slm))
final_vars <- names(which(coefs[-1,4] <= pval))

# Applicability Domain

















 






