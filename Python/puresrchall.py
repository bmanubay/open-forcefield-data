#-*- coding: utf-8 -*-
"""
Created on Wed May 18 11:35:57 2016

@author: bmanubay
"""
import thermopyl as th 
from thermopyl import thermoml_lib
import cirpy
import numpy as np
import pandas as pd
from sklearn.externals.joblib import Memory

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
   return cirpy.resolve(x, rtype)
   
# Compounds of most interest as decide by David, Chris and Bryce   
davmollist = ['2,2,4-trimethylpentane', 'cycloheptane', 'diisopropylether', 'isopropyl ether', 'dimethoxymethane', '2,3-dimethylbutane', '2,2-dimethylbutane', '3-methylpentane', 'neohexane', '4-methyl-2-pentanol', '2-methyl-2-pentanol', '1,1-diethoxyethane', 'tert-butanol', 'tetrahydrofuran', 'heptane', 'water', 'ethanol', '1-butanol', 'methyl tert-butyl ether']
S = pd.DataFrame({'IUPAC_Names': davmollist}, columns = ['IUPAC_Names'])
S["SMILES"] = S.IUPAC_Names.apply(lambda x: resolve_cached(x, "smiles"))

df = th.pandas_dataframe()
dt = list(df.columns)

bad_filenames = ["/home/bmanubay/.thermoml/j.fluid.2013.12.014.xml"]  # This file confirmed to have possible data entry errors.
df = df[~df.filename.isin(bad_filenames)]

experiments = ["Mass density, kg/m3","Speed of sound, m/s", "Relative permittivity at zero frequency", "Molar heat capacity at constant pressure, J/K/mol", "Molar enthalpy of vaporization or sublimation, kJ/mol", "Molar enthalpy, kJ/mol"]

ind_list = [df[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
df = df.ix[ind]

name_to_formula = pd.read_hdf("/home/bmanubay/.thermoml/compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()   


# Extract rows with two components
df["n_components"] = df.components.apply(lambda x: len(x.split("__")))
df = df[df.n_components == 1]
df.dropna(axis=1, how='all', inplace=True)

# Strip rows not in liquid phase
df = df[df['phase']=='Liquid']

df = df[df.components.isin(name_to_formula.index)]
name_to_formula = name_to_formula[name_to_formula.index.isin(df.components)]
df["formula"] = df.components.apply(lambda chemical: name_to_formula[chemical])

heavy_atoms = ["C", "O", "F", "N", "S", "B", "P", "Cl", "Br", "I"]
desired_atoms = ["H"] + heavy_atoms

df["n_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
df["n_heavy_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
df["n_desired_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
df["n_other_atoms"] = df.n_atoms - df.n_desired_atoms

df = df[df.n_other_atoms == 0]

df = df[df.n_heavy_atoms > 0]
df.dropna(axis=1, how='all', inplace=True)

df["SMILES"] = df.components.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.SMILES != None]
df.dropna(subset=["SMILES"], inplace=True)
df = df.ix[df.SMILES.dropna().index]

df["cas"] = df.components.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df["InChI"] = df.components.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "stdinchikey")))
df = df[df.cas != None]
df = df.ix[df.cas.dropna().index]

# Neither names (components) nor smiles are unique.  Use CAS to ensure consistency.
cannonical_smiles_lookup = df.groupby("cas").SMILES.first()
cannonical_components_lookup = df.groupby("cas").components.first()


df["SMILES"] = df.cas.apply(lambda x: cannonical_smiles_lookup[x])
df["components"] = df.cas.apply(lambda x: cannonical_components_lookup[x])

# Extract rows with temperature between 128 and 399 K
df = df[df['Temperature, K'] > 250.]
df = df[df['Temperature, K'] < 400.]

cols = ['filename', 'components', 'SMILES', 'cas', 'InChI', 'Temperature, K', 'Pressure, kPa', 'Molar enthalpy of vaporization or sublimation, kJ/mol', 'Molar enthalpy of vaporization or sublimation, kJ/mol_std']

dfhvap = df[cols]
dfhvap['Molar enthalpy of vaporization or sublimation, kJ/mol'].replace('nan', np.nan, inplace=True)
dfhvap = dfhvap[np.isnan(dfhvap['Molar enthalpy of vaporization or sublimation, kJ/mol'])==False]

# Extract rows with pressure between 101.325 kPa and 101325 kPa
df = df[df['Pressure, kPa'] > 100.]
df = df[df['Pressure, kPa'] < 102000.]

df.dropna(axis=1, how='all', inplace=True)

df["filename"] = df["filename"].map(lambda x: x.replace(' ', '')[25:])
df["filename"] = df.filename.map(lambda x: x.replace(' ', '')[:-4])

dfhvap["filename"] = dfhvap["filename"].map(lambda x: x.replace(' ', '')[25:])
dfhvap["filename"] = dfhvap.filename.map(lambda x: x.replace(' ', '')[:-4])

def dfpretty(df, prop):
    dfbig = pd.concat([df['filename'], df["components"], df["SMILES"], df["cas"], df["InChI"], df["Temperature, K"], df["Pressure, kPa"], df[prop], df[prop+"_std"]], axis=1, keys=["filename", "components", "SMILES", "CAS", "InChI", "Temperature, K", "Pressure, kPa", prop, prop+"_std"])
    dfbig[prop+"_std"].replace('nan', np.nan, inplace=True)
    dfbig = dfbig[np.isnan(dfbig[prop+"_std"])==False]
    cannonical_smiles_lookup = dfbig.groupby("CAS").SMILES.first()
    cannonical_components_lookup = dfbig.groupby("CAS").components.first()
    dfbig["SMILES"] = dfbig.CAS.apply(lambda x: cannonical_smiles_lookup[x])
    dfbig["components"] = dfbig.CAS.apply(lambda x: cannonical_components_lookup[x])
    a = dfbig["filename"].value_counts()
    a = a.reset_index()
    a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
    b = dfbig["InChI"].value_counts()
    b = b.reset_index()
    b.rename(columns={"index":"InChI","InChI":"Count"},inplace=True)
    b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name")) 
    b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles")) 
    
    return dfbig, a, b
    
                   
dfbig = pd.concat([df['filename'], df["components"], df["SMILES"], df["cas"], df["InChI"], df["Temperature, K"], df["Pressure, kPa"], df["Mass density, kg/m3"], df["Mass density, kg/m3_std"], df["Speed of sound, m/s"], df["Speed of sound, m/s_std"], df["Relative permittivity at zero frequency"], df["Relative permittivity at zero frequency_std"], df["Molar heat capacity at constant pressure, J/K/mol"], df["Molar heat capacity at constant pressure, J/K/mol_std"], df["Molar enthalpy, kJ/mol"], df["Molar enthalpy, kJ/mol_std"]] , axis=1, keys=["filename", "components", "SMILES", "CAS", "InChI", "Temperature, K", "Pressure, kPa", "Mass density, kg/m3", "Mass density, kg/m3_std", "Speed of sound, m/s", "Speed of sound, m/s_std", "Relative permittivity at zero frequency", "Relative permittivity at zero frequency_std", "Molar heat capacity at constant pressure, J/K/mol", "Molar heat capacity at constant pressure, J/K/mol_std", "Molar enthalpy, kJ/mol", "Molar enthalpy, kJ/mol_std"])
a = dfbig["filename"].value_counts()
a = a.reset_index()
a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b = dfbig["InChI"].value_counts()    
b = b.reset_index()    
b.rename(columns={"index":"InChI","InChI":"Count"},inplace=True)    
b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name"))     
b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles"))

df6 = dfhvap
a6 = df6["filename"].value_counts()
a6 = a6.reset_index()
a6.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b6 = df6["InChI"].value_counts()    
b6 = b6.reset_index()    
b6.rename(columns={"index":"InChI","InChI":"Count"},inplace=True)    
b6["Component"] = b6.InChI.apply(lambda x: resolve_cached(x, "iupac_name"))     
b6["SMILES"] = b6.InChI.apply(lambda x: resolve_cached(x, "smiles"))


df1, a1, b1 = dfpretty(df, "Mass density, kg/m3")
df2, a2, b2 = dfpretty(df, "Speed of sound, m/s")
df3, a3, b3 = dfpretty(df, "Relative permittivity at zero frequency")
df4, a4, b4 = dfpretty(df, "Molar heat capacity at constant pressure, J/K/mol")
df5, a5, b5 = dfpretty(df, "Molar enthalpy, kJ/mol")



pathdf = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Property data/"
pathjourn = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Journal name counts/"
pathcomp = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Component counts/"

def saveprettycsv(df, path, filename):
    df.to_csv(path+filename, sep =';')

def saveprettypickle(df, path, filename):
    df.to_pickle(path+filename)

# save csv with ; delimiter
saveprettycsv(dfbig, pathdf, "alldata_pure_full.csv")    
saveprettycsv(df1, pathdf, "dens_pure_full.csv") 
saveprettycsv(df2, pathdf, "sos_pure_full.csv") 
saveprettycsv(df3, pathdf, "dielec_pure_full.csv") 
saveprettycsv(df4, pathdf, "cpmol_pure_full.csv") 
saveprettycsv(df5, pathdf, "hmol_pure_full.csv") 
saveprettycsv(df6, pathdf, "hvap_pure_full.csv")

saveprettycsv(a, pathjourn, "purename_counts_all_full.csv")
saveprettycsv(a1, pathjourn, "purename_counts_dens_full.csv")
saveprettycsv(a2, pathjourn, "purename_counts_sos_full.csv")
saveprettycsv(a3, pathjourn, "purename_counts_dielec_full.csv")
saveprettycsv(a4, pathjourn, "purename_counts_cpmol_full.csv")
saveprettycsv(a5, pathjourn, "purename_counts_hmol_full.csv")
saveprettycsv(a6, pathjourn, "purename_counts_hvap_full.csv")

saveprettycsv(b, pathcomp, "purecomp_counts_all_full.csv")
saveprettycsv(b1, pathcomp, "purecomp_counts_dens_full.csv")
saveprettycsv(b2, pathcomp, "purecomp_counts_sos_full.csv")
saveprettycsv(b3, pathcomp, "purecomp_counts_dielec_full.csv")
saveprettycsv(b4, pathcomp, "purecomp_counts_cpmol_full.csv")
saveprettycsv(b5, pathcomp, "purecomp_counts_hmol_full.csv")
saveprettycsv(b6, pathcomp, "purecomp_counts_hvap_full.csv")

# save pickle
saveprettypickle(dfbig, pathdf, "alldata_pure_full.pkl")    
saveprettypickle(df1, pathdf, "dens_pure_full.pkl") 
saveprettypickle(df2, pathdf, "sos_pure_full.pkl") 
saveprettypickle(df3, pathdf, "dielec_pure_full.pkl") 
saveprettypickle(df4, pathdf, "cpmol_pure_full.pkl") 
saveprettypickle(df5, pathdf, "hmol_pure_full.pkl") 
saveprettypickle(df6, pathdf, "hvap_pure_full.pkl")

saveprettypickle(a, pathjourn, "purename_counts_all_full.pkl")
saveprettypickle(a1, pathjourn, "purename_counts_dens_full.pkl")
saveprettypickle(a2, pathjourn, "purename_counts_sos_full.pkl")
saveprettypickle(a3, pathjourn, "purename_counts_dielec_full.pkl")
saveprettypickle(a4, pathjourn, "purename_counts_cpmol_full.pkl")
saveprettypickle(a5, pathjourn, "purename_counts_hmol_full.pkl")
saveprettypickle(a6, pathjourn, "purename_counts_hvap_full.pkl")

saveprettypickle(b, pathcomp, "purecomp_counts_all_full.pkl")
saveprettypickle(b1, pathcomp, "purecomp_counts_dens_full.pkl")
saveprettypickle(b2, pathcomp, "purecomp_counts_sos_full.pkl")
saveprettypickle(b3, pathcomp, "purecomp_counts_dielec_full.pkl")
saveprettypickle(b4, pathcomp, "purecomp_counts_cpmol_full.pkl")
saveprettypickle(b5, pathcomp, "purecomp_counts_hmol_full.pkl") 
saveprettypickle(b6, pathcomp, "purecomp_counts_hvap_full.pkl")



