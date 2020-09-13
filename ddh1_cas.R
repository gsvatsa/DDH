library(rcdk, lib.loc = "lib")
library(itertools)

descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
moliter <- iload.molecules("antiviral_updated.sdf", type="sdf")
while(hasNext(moliter)) {
  mol <- nextElem(moliter)
  
  print(get.property(mol, "cdk:Title"))
}

mols <- sapply(smiles$SMILES, parse.smiles)
write.molecules(mols, "ddh01.sdf")

