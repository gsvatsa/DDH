#! /usr/bin/python

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

suppl = Chem.SDMolSupplier('/home/ChemTree/antiviral_updated.sdf')
df = pd.DataFrame()
for idx, mol in enumerate(suppl):
	try:
		name = mol.GetProp('cas.index.name')
		fr_ester = Descriptors.fr_ester(mol)
		fr_piperdine = Descriptors.fr_piperdine(mol)
		EState_VSA2 = Descriptors.EState_VSA2(mol)
		PEOE_VSA8 = Descriptors.PEOE_VSA8(mol)
	except:
		name = None
		fr_ester = None
		fr_piperdine = None
		EState_VSA2 = None
		PEOE_VSA8 = None
	finally:
		df2 = pd.DataFrame([{'name': name,
												'fr_ester': fr_ester,
                      	'fr_piperdine': fr_piperdine,
                      	'EState_VSA2': EState_VSA2,
                      	'PEOE_VSA8': PEOE_VSA8}], index = [idx + 1])
		print df2
		df = df.append(df2, ignore_index = True, sort = False)

df.to_csv("/home/ChemTree/cas_rdkit.tsv", sep = '\t')


