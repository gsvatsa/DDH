library(tidyverse)
library(rcdk)
library(matlib)

# Applicability Domain (AD)

all_rows = 1:nrow(smiles)
train_test_rows = which(!is.na(smiles$pIC50))
blinded_rows = all_rows[is.na(smiles$pIC50)]

# Hotelling's Test and Leverage (h)
# Pg.  37 (GUIDANCE DOCUMENT ON THE VALIDATION OF (QUANTITATIVE)STRUCTURE-ACTIVITY RELATIONSHIPS [(Q)SAR] MODELS)
library(Hotelling)
alpha = 0.05
train_moldesc <- moldesc[train_rows, model_vars]
test_moldesc <- moldesc[test_rows, model_vars]
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
desc <- train_moldesc
h_train <- apply(desc, 1, FUN = function(x) {
  t(x) %*% X_sq_inv %*% x
})s
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
  geom_text()

desc <- blinded_moldesc
h <- apply(desc, 1, FUN = function(x) {
  t(x) %*% X_sq_inv %*% x
})
# Warning Leverage
h_star <- 3*ncol(train_moldesc)/nrow(train_moldesc) 

in_AD_Leverage = which(h <= h_star)
out_AD_Leverage = which(h > h_star)

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

probs <- multivariate_chebychev(train_moldesc, blinded_moldesc)
in_AD_Chebychev <- blinded_rows[probs >= 0.05 & probs <= 0.95]
out_AD_Chebychev <- setdiff(blinded_rows, in_AD_rows)

# Tanimoto AD

train_mols <- parse.smiles(smiles$SMILES[train_test_rows])
test_mols <- parse.smiles(smiles$SMILES[blinded_rows])

fps_train <- lapply(train_mols, get.fingerprint, type='circular')
fps_test <- lapply(test_mols, get.fingerprint, type='circular')

fp.sim <- fingerprint::fp.sim.matrix(fps_train, fps_test, method='tanimoto')
max_rows <- apply(fp.sim, 2, which.max)










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

