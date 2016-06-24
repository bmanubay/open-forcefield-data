import numpy as np
import pandas as pd

# Read in property data produced from puresrch.py and binsrch.py
pathdfp = "/home/bmanubay/.thermoml/tables/Ken/Pure-Solvents/Property data/"
pathdfm = "/home/bmanubay/.thermoml/tables/Ken/Binary-Mixtures/Property data/"

a1 = pd.read_csv(pathdfp+"dens_pure.csv", sep=';')   
a2 = pd.read_csv(pathdfp+"hvap_pure.csv", sep=';')
a3 = pd.read_csv(pathdfp+"cpmol_pure.csv", sep=';')

b1 = pd.read_csv(pathdfm+"dens_bin.csv", sep=';')
b2 = pd.read_csv(pathdfm+"eme_bin.csv", sep=';')
b3 = pd.read_csv(pathdfm+"emcp_bin.csv", sep=';')

# Create new columns for combined SMILES strings on mixtures
#b1['SMILES_mix'] = b1['SMILES1'].map(str) + '__' + b1['SMILES2'].map(str)
#b2['SMILES_mix'] = b2['SMILES1'].map(str) + '__' + b2['SMILES2'].map(str)
#b3['SMILES_mix'] = b3['SMILES1'].map(str) + '__' + b3['SMILES2'].map(str)

# Merge property sets using filename as index
c1 = pd.merge(a1,b1,how='outer',on=['filename'])
c2 = pd.merge(a2,b2,how='outer',on=['filename'])
c3 = pd.merge(a3,b3,how='outer',on=['filename'])

print(c1["SMILES"])
