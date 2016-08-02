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
a0 = pd.read_csv(purcomppath+"purecomp_counts_all.csv", sep=';')
a1 = pd.read_csv(bincomppath+"bincomp_counts_all.csv", sep=';')
a2 = pd.read_csv(mixcomppath+"mix_counts_all.csv", sep=';')
a3 = pd.read_csv(purcomppath+"purecomp_counts_dens.csv", sep=';')
a4 = pd.read_csv(purcomppath+"purecomp_counts_sos.csv", sep=';')
a5 = pd.read_csv(purcomppath+"purecomp_counts_dielec.csv", sep=';')
a6 = pd.read_csv(purcomppath+"purecomp_counts_cpmol.csv", sep=';')
a7 = pd.read_csv(purcomppath+"purecomp_counts_hvap.csv", sep=';')
a8 = pd.read_csv(bincomppath+"bincomp_counts_dens.csv", sep=';')
a9 = pd.read_csv(bincomppath+"bincomp_counts_eme.csv", sep=';')
a10 = pd.read_csv(bincomppath+"bincomp_counts_emcp.csv", sep=';')
a11 = pd.read_csv(bincomppath+"bincomp_counts_emv.csv", sep=';')
a12 = pd.read_csv(bincomppath+"bincomp_counts_actcoeff.csv", sep=';')
a13 = pd.read_csv(bincomppath+"bincomp_counts_sos.csv", sep=';')
a14 = pd.read_csv(bincomppath+"bincomp_counts_dielec.csv", sep=';')
a15 = pd.read_csv(mixcomppath+"mix_counts_dens.csv", sep=';')
a16 = pd.read_csv(mixcomppath+"mix_counts_eme.csv", sep=';')
a17 = pd.read_csv(mixcomppath+"mix_counts_emcp.csv", sep=';')
a18 = pd.read_csv(mixcomppath+"mix_counts_emv.csv", sep=';')
a19 = pd.read_csv(mixcomppath+"mix_counts_actcoeff.csv", sep=';')
a20 = pd.read_csv(mixcomppath+"mix_counts_sos.csv", sep=';')
a21 = pd.read_csv(mixcomppath+"mix_counts_dielec.csv", sep=';')

b = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_chain_filt1.smi",delim_whitespace=True,header=None)
c = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_rings_filt1.smi",delim_whitespace=True,header=None)

def mixsplit(df):
    df['x1'], df['x2'] = zip(*df['Mixture'].str.split('__').tolist())
    df['x2'].replace('', np.nan, inplace=True)
    df.dropna(subset=['x2'], inplace=True) 

    df['SMILES1'] = df.x1.apply(lambda x: resolve_cached(x, 'smiles'))
    df['SMILES2'] = df.x2.apply(lambda x: resolve_cached(x, 'smiles'))

    return df

a2 = mixsplit(a2)
a15 = mixsplit(a15)
a16 = mixsplit(a16)
a17 = mixsplit(a17)
a18 = mixsplit(a18)
a19 = mixsplit(a19)
a20 = mixsplit(a20)
a21 = mixsplit(a21)


d = pd.merge(b,c,on=0,how='outer')
d.columns = ["SMILES","chains,","rings"]

# Return index list of booleans for which components of read count df is in Chris's list
ind0 = a0.SMILES.isin(d.SMILES)
ind1 = a1.SMILES.isin(d.SMILES)
ind2 = a2.SMILES1.isin(d.SMILES) & a2.SMILES2.isin(d.SMILES) 
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
ind15 = a15.SMILES1.isin(d.SMILES) & a15.SMILES2.isin(d.SMILES)
ind16 = a16.SMILES1.isin(d.SMILES) & a16.SMILES2.isin(d.SMILES)
ind17 = a17.SMILES1.isin(d.SMILES) & a17.SMILES2.isin(d.SMILES)
ind18 = a18.SMILES1.isin(d.SMILES) & a18.SMILES2.isin(d.SMILES)
ind19 = a19.SMILES1.isin(d.SMILES) & a19.SMILES2.isin(d.SMILES)
ind20 = a20.SMILES1.isin(d.SMILES) & a20.SMILES2.isin(d.SMILES)
ind21 = a21.SMILES1.isin(d.SMILES) & a21.SMILES2.isin(d.SMILES)

# Create new df based on boolean lists
a0_diverse_chk = a0[ind0]
a0_diverse_chk = a0_diverse_chk.drop("Unnamed: 0", axis = 1)
a1_diverse_chk = a1[ind1]
a1_diverse_chk = a1_diverse_chk.drop("Unnamed: 0", axis = 1)
a2_diverse_chk = a2[ind2]
a2_diverse_chk = a2_diverse_chk.drop("Unnamed: 0", axis = 1)
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
a15_diverse_chk = a15[ind15]
a15_diverse_chk = a15_diverse_chk.drop("Unnamed: 0", axis = 1)
a16_diverse_chk = a16[ind16]
a16_diverse_chk = a16_diverse_chk.drop("Unnamed: 0", axis = 1)
a17_diverse_chk = a17[ind17]
a17_diverse_chk = a17_diverse_chk.drop("Unnamed: 0", axis = 1)
a18_diverse_chk = a18[ind18]
a18_diverse_chk = a18_diverse_chk.drop("Unnamed: 0", axis = 1)
a19_diverse_chk = a19[ind19]
a19_diverse_chk = a19_diverse_chk.drop("Unnamed: 0", axis = 1)
a20_diverse_chk = a20[ind20]
a20_diverse_chk = a20_diverse_chk.drop("Unnamed: 0", axis = 1)
a21_diverse_chk = a21[ind21]
a21_diverse_chk = a21_diverse_chk.drop("Unnamed: 0", axis = 1)

# save diverse checked df as csv and pickle files
a0_diverse_chk.to_csv(purcomppath+"purecomp_counts_all_AlkEthOH.csv", sep=';') 
a1_diverse_chk.to_csv(bincomppath+"bincomp_counts_all_AlkEthOH.csv", sep=';')
a2_diverse_chk.to_csv(mixcomppath+"mix_counts_all_AlkEthOH.csv", sep=';')
a3_diverse_chk.to_csv(purcomppath+"purecomp_counts_dens_AlkEthOH.csv", sep=';')
a4_diverse_chk.to_csv(purcomppath+"purecomp_counts_sos_AlkEthOH.csv", sep=';')
a5_diverse_chk.to_csv(purcomppath+"purecomp_counts_dielec_AlkEthOH.csv", sep=';')
a6_diverse_chk.to_csv(purcomppath+"purecomp_counts_cpmol_AlkEthOH.csv", sep=';')
a7_diverse_chk.to_csv(purcomppath+"purecomp_counts_hvap_AlkEthOH.csv", sep=';')
a8_diverse_chk.to_csv(bincomppath+"bincomp_counts_dens_AlkEthOH.csv", sep=';')
a9_diverse_chk.to_csv(bincomppath+"bincomp_counts_eme_AlkEthOH.csv", sep=';')
a10_diverse_chk.to_csv(bincomppath+"bincomp_counts_emcp_AlkEthOH.csv", sep=';')
a11_diverse_chk.to_csv(bincomppath+"bincomp_counts_emv_AlkEthOH.csv", sep=';')
a12_diverse_chk.to_csv(bincomppath+"bincomp_counts_actcoeff_AlkEthOH.csv", sep=';')
a13_diverse_chk.to_csv(bincomppath+"bincomp_counts_sos_AlkEthOH.csv", sep=';')
a14_diverse_chk.to_csv(bincomppath+"bincomp_counts_dielec_AlkEthOH.csv", sep=';')
a15_diverse_chk.to_csv(mixcomppath+"mix_counts_dens_AlkEthOH.csv", sep=';')
a16_diverse_chk.to_csv(mixcomppath+"mix_counts_eme_AlkEthOH.csv", sep=';')
a17_diverse_chk.to_csv(mixcomppath+"mix_counts_emcp_AlkEthOH.csv", sep=';')
a18_diverse_chk.to_csv(mixcomppath+"mix_counts_emv_AlkEthOH.csv", sep=';')
a19_diverse_chk.to_csv(mixcomppath+"mix_counts_actcoeff_AlkEthOH.csv", sep=';')
a20_diverse_chk.to_csv(mixcomppath+"mix_counts_sos_AlkEthOH.csv", sep=';')
a21_diverse_chk.to_csv(mixcomppath+"mix_counts_dielec_AlkEthOH.csv", sep=';')

a0_diverse_chk.to_pickle(purcomppath+"purecomp_counts_all_AlkEthOH.pkl") 
a1_diverse_chk.to_pickle(bincomppath+"bincomp_counts_all_AlkEthOH.pkl")
a2_diverse_chk.to_pickle(mixcomppath+"mix_counts_all_AlkEthOH.pkl")

