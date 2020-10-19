library(rcdk, lib.loc = "lib")
library(iterators)
library(tidyverse)

cas_smiles_file <- "results/cas_smiles.csv"
cas_sdf_file <- "antiviral_updated.sdf"
if (!file.exists(cas_smiles_file)) {
  moliter <- iload.molecules(cas_sdf_file, type="sdf")
  cas_smiles <- character()
  while(rcdk::hasNext(moliter)) {
    mol <- nextElem(moliter)
    cas_smiles <- append(cas_smiles, get.smiles(mol))
  }
  cas_smiles_df <-
    data.frame(SMILES = cas_smiles) %>%
    mutate(CAS_ID = row_number()) %>%
    relocate(CAS_ID) %>%
    write_excel_csv(cas_smiles_file)
} else {
  cas_smiles_df <- read_csv(cas_smiles_file, 
                            col_types = cols(CAS_ID = col_integer()))
}

padel_model <- models[["padel"]]
padel_model_vars <- names(coef(padel_model))[-1]
cas_padel <- 
  read_csv("cas_padel2.zip") %>% 
  mutate(CAS_ID = row_number()) %>%
  filter_all(all_vars(!is.infinite(.) & !is.na(.))) #%>%
  #predict(scalers[["padel"]], newdata = .)

df1 <- data.frame(CAS_ID = cas_padel$CAS_ID) 
df1$padel = predict(padel_model, 
                    newdata = predict(scalers[["padel"]], cas_padel))

cas_rdkit <- 
  read_tsv("cas_rdkit.tsv") %>% 
  mutate(CAS_ID = row_number()) %>%
  filter_all(all_vars(!is.infinite(.) & !is.na(.)))

cas_padel_rdkit <- inner_join(cas_padel, cas_rdkit, by="CAS_ID") 
df2 <- data.frame(CAS_ID = cas_padel_rdkit$CAS_ID)
df2$padel_rdkit <- predict(models[["padel_rdkit"]], 
                           newdata = predict(scalers[["padel_rdkit"]], cas_padel_rdkit))

# CAS Predictions
# preds_cas_ensemble <- predict(rf_fit, 
#                               data.frame(padel = preds_cas_padel, 
#                                          padel_rdkit = preds_cas_padel_rdkit))
cas_df <- 
  inner_join(df1, df2, by = "CAS_ID") %>%
  filter_all(all_vars(!is.infinite(.) & !is.na(.))) %>%
  mutate(pIC50_pred = predict(rf_fit, newdata = .))
  

# CAS AD
cas_moldescs = list()
cas_moldescs[["padel"]] = cas_padel#[cas_df$CAS_ID, ]
cas_moldescs[["padel_rdkit"]] = cas_padel_rdkit#[cas_df$CAS_ID,]

cas_AD_Infos <- lapply(sel_desc_types, FUN = function(desc_type) {
  moldesc <- read_csv(paste0("descriptors/", desc_type, "_desc.csv")) 
  model <- SMLR_models[[desc_type]]
  model_vars <- gsub('\`', '', variable.names(model)[-1])
  scaler <- scalers[[desc_type]]
  
  train_moldesc <- 
    moldesc[train_rows,] %>% 
    predict(scaler, newdata = .) %>%
    select(all_of(model_vars))
  cas_moldesc <- 
    cas_moldescs[[desc_type]] %>%
    predict(scaler, newdata = .) %>%
    select(all_of(model_vars))
  
  AD_Leverage_cas <- Leverage_AD(model, train_moldesc, cas_moldesc, 
                                 type = "cas", plot_train = F)
  AD_Chebychev_cas <- Chebychev_AD(train_moldesc, cas_moldesc)
 
  AD_Standardization_cas <- Standardization_AD(desc_type,
                                               train_moldesc, 
                                               cas_moldesc,
                                               type = "cas")
  
  AD_Extrapolation_cas <- Extrapolation_AD(desc_type, cas_moldesc)
                                           
  
  return(list(AD_Leverage_cas = AD_Leverage_cas,
              AD_Chebychev_cas = AD_Chebychev_cas,
              AD_Standardization_cas = AD_Standardization_cas,
              AD_Extrapolation_cas = AD_Extrapolation_cas))
})
names(cas_AD_Infos) <- sel_desc_types


out_AD_cas <- Reduce("intersect", lapply(sel_desc_types, FUN = function(desc_type) {
  AD <- cas_AD_Infos[[desc_type]]
  out_AD_Leverage_cas <- AD[["AD_Leverage_cas"]][[2]]
  out_AD_Chebychev_cas <- AD[["AD_Chebychev_cas"]][[2]]
  out_AD_Standardization_cas <- AD[["AD_Standardization_cas"]][[2]]
  out_AD_Extrapolation_cas <- AD[["AD_Extrapolation_cas"]][[2]]
  
  cas_ids <- cas_moldescs[[desc_type]]$CAS_ID
  Reduce("union", list(cas_ids[out_AD_Leverage_cas],
                       cas_ids[out_AD_Chebychev_cas],
                       cas_ids[out_AD_Standardization_cas],
                       cas_ids[out_AD_Extrapolation_cas]
  ))
}))


# DDT Input Form 3
cas_df_100 <- cas_df %>%
  arrange(desc(pIC50_pred)) %>%
  head(100) %>%
  mutate(SMILES = cas_smiles_df$SMILES[CAS_ID]) %>%
  relocate(.after = CAS_ID) %>%
  mutate(AD = if_else(CAS_ID %in% out_AD_cas, "Outside_AD", "-")) %>%
  select(CAS_ID, SMILES, pIC50_pred, AD) %>%
  write_excel_csv("results/Input Form 3_DDT1-01.csv")




