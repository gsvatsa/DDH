# Drug Discovery Hackathon - Developing a regression-based model for screening compounds with 3CLpro inhibitory activity


A data set of diverse compounds with 3CLpro inhibitory activity (pIC50 with IC50 in microM) is provided. The SMILES notations of all compounds are given. The activity of a few compounds is blinded. The remaining set of compounds should be divided into a training set and a test set using 70:30 ratios. Simple, reproducible and easily transferable regression-based QSAR models should be developed from the training set compounds using only 2D descriptors. The models should be validated based on the test set compounds. The models should have at least the following quality: R2 > 0.7, LOO-Q2 > 0.7, Q2ext_F1 > 0.7, Tropsha’s criteria: Pass, MAE based criteria: Moderate or good. Using the best model, the activity of the blinded set should be determined. Care should also be taken to consider the applicability domain while predicting the external/blinded compounds. The best model should also be used to prioritize CAS antiviral database compounds for 3CLpro inhibitory activity and the top ranking 100 compounds should be listed (along with AD information). Additional Comments: i) The participants may explore new methods of feature/descriptor selection, new methods of selection of appropriate training set compounds, new methods of defining applicability domain. ii) The participants may collect from the literature new data with the experimental activity/property values and show the applicability of the developed model, and iii) The participants may explore the response property or any other related measure computed from another method for the prioritized compounds and show consistency with the predicted values from the developed model. A consensus will further prove reliability of the developed model.

# Code

Linear Regression model development and QSAR validation code is in ddh1_lm.R

Application to CAS database and final results in ddh1_cas.R

