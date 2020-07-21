library(tidyverse)
library(rcdk)
library(matlib)

# Applicability Domain (AD)

all_rows = 1:nrow(smiles)
train_test_rows = which(!is.na(smiles$pIC50))
blinded_rows = all_rows[is.na(smiles$pIC50)]

set.seed(1024)
test_rows = sample(train_test_rows, 0.3*floor(length(train_test_rows)))
train_rows = setdiff(train_test_rows, test_rows)

# Hotelling's Test and Leverage (h)
# Pg.  37 (GUIDANCE DOCUMENT ON THE VALIDATION OF (QUANTITATIVE)STRUCTURE-ACTIVITY RELATIONSHIPS [(Q)SAR] MODELS)
library(Hotelling)
alpha = 0.05
train_moldesc <- moldesc[train_rows, model_vars]
test_moldesc <- moldesc[test_rows, model_vars]
train_test_moldesc <- moldesc[train_test_rows, model_vars]
blinded_moldesc <- moldesc[blinded_rows, model_vars]

fit <- hotelling.test(train_test_moldesc, blinded_moldesc)
if (fit$pval > alpha) {
  print("H0: Blinded and Train mol samples have the same mean")
} else {
  print("H1: Blinded and Train mol samples have different means")
}

# Leverage and Williams Plot
X <- as.matrix(train_moldesc)
X_sq_inv <- Inverse(t(X) %*% X)

# Warning Leverage
h_warn <- 3*(ncol(X) + 1)/nrow(X) 

# Williams Plot
desc <- train_moldesc
h_train <- apply(desc, 1, FUN = function(x) {
  t(x) %*% X_sq_inv %*% x
})
w_train <- data.frame(stdres = train_stdres, h = h_train, class = "train")
rownames(w_train) <- paste("T",seq_along(h_train), sep="")
desc <- test_moldesc
h_test <- apply(desc, 1, FUN = function(x) {
  t(x) %*% X_sq_inv %*% x
})
w_test <- data.frame(stdres = test_stdres, h = h_test, class = "test")
rownames(w_test) <- paste("S",seq_along(h_test), sep="")
w <- rbind(w_train, w_test)

ggplot(data = w, aes(x = h, y = stdres, shape = class, label = rownames(w))) +
  geom_point() +
  geom_text() +
  geom_hline(yintercept = mean(w_train$stdres) + 2*sd(w_train$stdres), color = "red") +
  geom_hline(yintercept = mean(w_train$stdres) - 2*sd(w_train$stdres), color = "red") +
  geom_vline(xintercept = h_warn, color = "red")

desc <- blinded_moldesc
h <- apply(desc, 1, FUN = function(x) {
  t(x) %*% X_sq_inv %*% x
})

in_AD_test_Leverage = which(h_test <= h_warn)
out_AD_test_Leverage = which(h_test > h_warn)

in_AD_Leverage = which(h <= h_warn)
out_AD_Leverage = which(h > h_warn)

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

train_test_moldesc <- moldesc[train_test_rows, model_vars]
train_moldesc <- moldesc[train_rows, model_vars]
test_moldesc <- moldesc[test_rows, model_vars]
blinded_moldesc <- moldesc[blinded_rows, model_vars]

probs_train <- multivariate_chebychev(train_moldesc, train_moldesc)
probs_test <- multivariate_chebychev(train_moldesc, test_moldesc)
probs_train_test <- multivariate_chebychev(train_test_moldesc, train_test_moldesc)
probs <- multivariate_chebychev(train_moldesc, blinded_moldesc)

in_AD_test_Chebychev <- which(probs_test >= min(probs_train) & probs_test <= max(probs_train))
out_AD_test_Chebychev <- seq_along(probs_test)[-in_AD_Chebychev]

in_AD_Chebychev <- which(probs >= min(probs_train) & probs <= max(probs_train))
out_AD_Chebychev <- seq_along(probs)[-in_AD_Chebychev]

# Tanimoto AD

train_mols <- parse.smiles(smiles$SMILES[train_rows])
test_mols <- parse.smiles(smiles$SMILES[test_rows])
train_test_mols <- parse.smiles(smiles$SMILES[train_test_rows])
blinded_mols <- parse.smiles(smiles$SMILES[blinded_rows])

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
sig_types = c('standard', 'extended', 'graph', 'hybridization', 
              'maccs', 'estate', 'pubchem', 'shortestpath', 
              'circular', 'substructure')
options("java.parameters"=c("-Xmx4000m"))
rmse <- sapply(sig_types, FUN = function(sig_type) {
  
  fps_train <- lapply(train_mols, get.fingerprint, type=sig_type)
  fps_test <- lapply(test_mols, get.fingerprint, type=sig_type)
  test_sim <- fingerprint::fp.sim.matrix(fps_train, fps_test, method='tanimoto')
  
  # 1-NN Tanimoto
  tanimoto_1nn <- apply(test_sim, 2, which.max)
  tanimoto_1nn_sim <- sapply(seq_along(tanimoto_1nn), FUN = function(i) {
    y <- y_test[i]
    row <- tanimoto_1nn[i]
    abs(y_train[row] - y)
  })

  sqrt(mean((activ_1nn_sim - tanimoto_1nn_sim)^2))

})

sig_type = names(which.min(rmse))
fps_train <- lapply(train_mols, get.fingerprint, type=sig_type)
fps_blinded <- lapply(blinded_mols, get.fingerprint, type=sig_type)

blinded_sim <- fingerprint::fp.sim.matrix(fps_train, fps_blinded, method='tanimoto')
tanimoto_1nn <- apply(blinded_sim, 2, which.max)
tanimoto_1nn_sim <- sapply(seq_along(tanimoto_1nn), FUN = function(col) {
  blinded_sim[tanimoto_1nn[col], col]
})

sim_threshold <- 0.8

in_AD_Tanimoto <- which(tanimoto_1nn_sim >= sim_threshold)
out_AD_Tanimoto <- which(tanimoto_1nn_sim < sim_threshold)

out_AD <- Reduce(intersect, list(out_AD_Leverage, out_AD_Chebychev, out_AD_Tanimoto))















library(readxl)
test_AD <- read_excel("~/Downloads/AD-MDI version 1.2/Output/blinded_AD.xlsx")
names(test_AD)[3] = "Pred_Error"

fit <- lm(Pred_Error ~ MDI - 1, data = test_AD)
summary(fit)
plot(Pred_Error ~ MDI, data = test_AD)

library(readr)
blinded_QueryMDI <- read_delim("~/Downloads/AD-MDI version 1.2/Output/blinded_QueryMDI.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE, 
                               skip = 1)
blinded_QueryMDI$Pred_Error <- predict(fit, newdata = blinded_QueryMDI)

