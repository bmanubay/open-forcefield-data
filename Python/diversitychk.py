# -*- coding: utf-8 -*-
"""
Created on Sun May  8 18:29:53 2016

@author: bmanubay
"""

# Check what moelcules we have appear in Chris's list

import pandas as pd
import numpy as np
import cirpy
from sklearn.externals.joblib import Memory

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
    return cirpy.resolve(x, rtype)

purcomppath = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Component counts/"
bincomppath = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Binary-Mixtures/Component counts/"
mixcomppath = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Binary-Mixtures/Mixture counts/"

# read in ; delimited csv of comp/mix counts created in thermomlcnts.py
a0 = pd.read_csv(purcomppath+"purecomp_counts_all_full.csv", sep=';')
a1 = pd.read_csv(bincomppath+"bincomp_counts_all_full.csv", sep=';')
a3 = pd.read_csv(purcomppath+"purecomp_counts_dens_full.csv", sep=';')
a4 = pd.read_csv(purcomppath+"purecomp_counts_sos_full.csv", sep=';')
a5 = pd.read_csv(purcomppath+"purecomp_counts_dielec_full.csv", sep=';')
a6 = pd.read_csv(purcomppath+"purecomp_counts_cpmol_full.csv", sep=';')
a7 = pd.read_csv(purcomppath+"purecomp_counts_hvap_full.csv", sep=';')
a8 = pd.read_csv(bincomppath+"bincomp_counts_dens_full.csv", sep=';')
a9 = pd.read_csv(bincomppath+"bincomp_counts_eme_full.csv", sep=';')
a10 = pd.read_csv(bincomppath+"bincomp_counts_emcp_full.csv", sep=';')
a11 = pd.read_csv(bincomppath+"bincomp_counts_emv_full.csv", sep=';')
a12 = pd.read_csv(bincomppath+"bincomp_counts_actcoeff_full.csv", sep=';')
a13 = pd.read_csv(bincomppath+"bincomp_counts_sos_full.csv", sep=';')
a14 = pd.read_csv(bincomppath+"bincomp_counts_dielec_full.csv", sep=';')

b = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_chain_filt1.smi",delim_whitespace=True,header=None)
c = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_rings_filt1.smi",delim_whitespace=True,header=None)



d = pd.merge(b,c,on=0,how='outer')
d.columns = ["SMILES","chains,","rings"]

# Return index list of booleans for which components of read count df is in Chris's list
ind0 = a0.SMILES.isin(d.SMILES)
ind1 = a1.SMILES.isin(d.SMILES)

ind3 = a3.SMILES.isin(d.SMILES)
ind4 = a4.SMILES.isin(d.SMILES)
ind5 = a5.SMILES.isin(d.SMILES)
ind6 = a6.SMILES.isin(d.SMILES)
ind7 = a7.SMILES.isin(d.SMILES)
ind8 = a8.SMILES.isin(d.SMILES)
ind9 = a9.SMILES.isin(d.SMILES)
ind10 = a10.SMILES.isin(d.SMILES)
ind11 = a11.SMILES.isin(d.SMILES)
ind12 = a12.SMILES.isin(d.SMILES)
ind13 = a13.SMILES.isin(d.SMILES)
ind14 = a14.SMILES.isin(d.SMILES)

# Create new df based on boolean lists
a0_diverse_chk = a0[ind0]
a0_diverse_chk = a0_diverse_chk.drop("Unnamed: 0", axis = 1)
a1_diverse_chk = a1[ind1]
a1_diverse_chk = a1_diverse_chk.drop("Unnamed: 0", axis = 1)
a3_diverse_chk = a3[ind3]
a3_diverse_chk = a3_diverse_chk.drop("Unnamed: 0", axis = 1)
a4_diverse_chk = a4[ind4]
a4_diverse_chk = a4_diverse_chk.drop("Unnamed: 0", axis = 1)
a5_diverse_chk = a5[ind5]
a5_diverse_chk = a5_diverse_chk.drop("Unnamed: 0", axis = 1)
a6_diverse_chk = a6[ind6]
a6_diverse_chk = a6_diverse_chk.drop("Unnamed: 0", axis = 1)
a7_diverse_chk = a7[ind7]
a7_diverse_chk = a7_diverse_chk.drop("Unnamed: 0", axis = 1)
a8_diverse_chk = a8[ind8]
a8_diverse_chk = a8_diverse_chk.drop("Unnamed: 0", axis = 1)
a9_diverse_chk = a9[ind9]
a9_diverse_chk = a9_diverse_chk.drop("Unnamed: 0", axis = 1)
a10_diverse_chk = a10[ind10]
a10_diverse_chk = a10_diverse_chk.drop("Unnamed: 0", axis = 1)
a11_diverse_chk = a11[ind11]
a11_diverse_chk = a11_diverse_chk.drop("Unnamed: 0", axis = 1)
a12_diverse_chk = a12[ind12]
a12_diverse_chk = a12_diverse_chk.drop("Unnamed: 0", axis = 1)
a13_diverse_chk = a13[ind13]
a13_diverse_chk = a13_diverse_chk.drop("Unnamed: 0", axis = 1)
a14_diverse_chk = a14[ind14]
a14_diverse_chk = a14_diverse_chk.drop("Unnamed: 0", axis = 1)

# save diverse checked df as csv and pickle files
a0_diverse_chk.to_csv(purcomppath+"purecomp_counts_all_full_AlkEthOH.csv", sep=';') 
a1_diverse_chk.to_csv(bincomppath+"bincomp_counts_all_full_AlkEthOH.csv", sep=';')
a3_diverse_chk.to_csv(purcomppath+"purecomp_counts_dens_full_AlkEthOH.csv", sep=';')
a4_diverse_chk.to_csv(purcomppath+"purecomp_counts_sos_full_AlkEthOH.csv", sep=';')
a5_diverse_chk.to_csv(purcomppath+"purecomp_counts_dielec_full_AlkEthOH.csv", sep=';')
a6_diverse_chk.to_csv(purcomppath+"purecomp_counts_cpmol_full_AlkEthOH.csv", sep=';')
a7_diverse_chk.to_csv(purcomppath+"purecomp_counts_hvap_full_AlkEthOH.csv", sep=';')
a8_diverse_chk.to_csv(bincomppath+"bincomp_counts_dens_full_AlkEthOH.csv", sep=';')
a9_diverse_chk.to_csv(bincomppath+"bincomp_counts_eme_full_AlkEthOH.csv", sep=';')
a10_diverse_chk.to_csv(bincomppath+"bincomp_counts_emcp_full_AlkEthOH.csv", sep=';')
a11_diverse_chk.to_csv(bincomppath+"bincomp_counts_emv_full_AlkEthOH.csv", sep=';')
a12_diverse_chk.to_csv(bincomppath+"bincomp_counts_actcoeff_full_AlkEthOH.csv", sep=';')
a13_diverse_chk.to_csv(bincomppath+"bincomp_counts_sos_full_AlkEthOH.csv", sep=';')
a14_diverse_chk.to_csv(bincomppath+"bincomp_counts_dielec_full_AlkEthOH.csv", sep=';')

a0_diverse_chk.to_pickle(purcomppath+"purecomp_counts_all_full_AlkEthOH.pkl") 
a1_diverse_chk.to_pickle(bincomppath+"bincomp_counts_all_full_AlkEthOH.pkl")

