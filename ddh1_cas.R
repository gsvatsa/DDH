library(rcdk, lib.loc = "lib")
library(itertools)

descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
moliter <- iload.molecules("antiviral_updated.sdf", type="sdf")
cas_smiles <- character()
while(hasNext(moliter)) {
  mol <- nextElem(moliter)
  cas_smiles <- append(cas_smiles, get.smiles(mol))
}

mols <- sapply(smiles$SMILES, parse.smiles)
write.molecules(mols, "ddh01.sdf")

cas_padel <- read_csv("cas_padel2.zip")
preds_cas_padel <- predict(models[["padel"]], newdata = cas_padel)

cas_rdkit <- read_tsv("cas_rdkit.tsv")

cas_padel_rdkit <- cbind(cas_padel, cas_rdkit)
preds_cas_padel_rdkit <- predict(models[["padel_rdkit"]], newdata = cas_padel_rdkit)

# CAS Predictions
preds_cas_ensemble <- predict(rf_fit, 
                              data.frame(padel = preds_cas_padel, 
                                         padel_rdkit = preds_cas_padel_rdkit))
cas_df <- 
  data.frame(padel = preds_cas_padel, padel_rdkit = preds_cas_padel_rdkit) %>%
  mutate(CAS_ID = 1:nrow(.)) %>%
  relocate(CAS_ID) %>%
  filter(!(is.na(preds_cas_padel) | is.na(preds_cas_padel_rdkit))) %>%
  mutate(pIC50_pred = predict(rf_fit, newdata = .)) %>%
  arrange(desc(pIC50_pred)) %>%
  head(100) %>%
  mutate(SMILES = cas_smiles[CAS_ID]) %>%
  relocate(.after = CAS_ID)

# CAS AD
cas_moldescs = list()
cas_moldescs[["padel"]] = cas_padel[cas_df$CAS_ID, ]
cas_moldescs[["padel_rdkit"]] = cas_padel_rdkit[cas_df$CAS_ID,]

cas_AD_Infos <- lapply(sel_desc_types, FUN = function(desc_type) {
  moldesc <- read_csv(paste0("descriptors/", desc_type, "_desc.csv")) 
  model <- SMLR_models[[desc_type]]
  model_vars <- gsub('\`', '', variable.names(model)[-1])
  
  train_moldesc <- moldesc[train_rows, model_vars]
  train_test_moldesc <- moldesc[train_test_rows, model_vars]
  cas_moldesc <- cas_moldescs[[desc_type]][, model_vars]
  
  AD_Leverage_cas <- Leverage_AD(model, train_moldesc, cas_moldesc, 
                                 type = "cas", plot_train = F)
  AD_Chebychev_cas <- Chebychev_AD(train_moldesc, cas_moldesc)
  # TODO 
  # AD_Standardization_cas <- Standardization_AD(train_test_moldesc, 
  #                                              cas_moldesc,
  #                                              type = "cas")
  return(list(AD_Leverage_cas = AD_Leverage_cas,
              AD_Chebychev_cas = AD_Chebychev_cas))
})
names(cas_AD_Infos) <- sel_desc_types

out_AD_cas <- Reduce("intersect", lapply(sel_desc_types, FUN = function(desc_type) {
  AD <- cas_AD_Infos[[desc_type]]
  out_AD_Leverage_cas <- AD[["AD_Leverage_cas"]][[2]]
  out_AD_Chebychev_cas <- AD[["AD_Chebychev_cas"]][[2]]
  #out_AD_Standardization_cas <- AD[["AD_Standardization_cas"]][[2]]
  
  Reduce("union", list(out_AD_Leverage_cas,
                       out_AD_Chebychev_cas
                       #out_AD_Standardization_blinded
  ))
}))

# DDT Input Form 3
cas_df %>%
  mutate(AD = ifelse(row_number() %in% out_AD_cas, "Outside_AD", "-")) %>%
  select(CAS_ID, SMILES, pIC50_pred, AD) %>%
  write_excel_csv("results/Input Form 3_DDT1-01.csv")




