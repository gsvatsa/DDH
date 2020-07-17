library(tidyverse)
library(rcdk)
library(matlib)

# Applicability Domain (AD)

train_test_mols <- moldesc[train_test_rows, model_vars]
blinded_mols <- moldesc[blinded_rows, model_vars]

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

probs <- multivariate_chebychev(train_test_mols, blinded_mols)
probs_test <- multivariate_chebychev(train_test_mols, train_test_mols)
in_AD <- blinded_rows[probs >= 0.05 & probs <= 0.95]
out_AD <- setdiff(blinded_rows, in_AD_rows)








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

