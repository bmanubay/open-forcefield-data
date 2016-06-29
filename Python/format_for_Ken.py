import numpy as np
import pandas as pd

# Read in property data produced from puresrch.py and binsrch.py
pathdfp = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Property data/"
pathdfm = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Binary-Mixtures/Property data/"

a1 = pd.read_csv(pathdfp+"dens_pure.csv", sep=';', usecols=['filename', 'components'])   
a2 = pd.read_csv(pathdfp+"hvap_pure.csv", sep=';', usecols=['filename', 'components'])
a3 = pd.read_csv(pathdfp+"cpmol_pure.csv", sep=';', usecols=['filename', 'components'])

b1 = pd.read_csv(pathdfm+"dens_bin.csv", sep=';', usecols=['filename', 'components'])
b2 = pd.read_csv(pathdfm+"eme_bin.csv", sep=';', usecols=['filename', 'components'])
b3 = pd.read_csv(pathdfm+"emcp_bin.csv", sep=';', usecols=['filename', 'components'])

# Create new columns for combined SMILES strings on mixtures
b1 = b1.rename(columns = {'components':'components_mix'})
b2 = b2.rename(columns = {'components':'components_mix'})
b3 = b3.rename(columns = {'components':'components_mix'})

# Find unique compoounds/mixtures for each journal in datasets 
d1 = a1.groupby('filename').components.unique()
d1 = d1.reset_index()
d1.columns = ['filename', 'Mass density, kg/m3']
e1 = b1.groupby('filename').components_mix.unique()
e1 = e1.reset_index()
e1.columns = ['filename', 'Mass density, kg/m3_binary']
d2 = a2.groupby('filename').components.unique()
d2 = d2.reset_index()
d2.columns = ['filename', 'Molar enthalpy of vaporization or sublimation, kJ/mol']
e2 = b2.groupby('filename').components_mix.unique()
e2 = e2.reset_index()
e2.columns = ['filename', 'Excess molar enthalpy (molar enthalpy of mixing), kJ/mol']
d3 = a3.groupby('filename').components.unique()
d3 = d3.reset_index()
d3.columns = ['filename', 'Molar heat capacity at constant pressure, J/K/mol']
e3 = b3.groupby('filename').components_mix.unique()
e3 = e3.reset_index()
e3.columns = ['filename', 'Excess molar heat capacity, J/K/mol']

# Merge the unique component/mixture list with filename as index
f1 = pd.merge(d1,e1,how='outer',on=['filename'])
f2 = pd.merge(d2,e2,how='outer',on=['filename'])
f3 = pd.merge(d3,e3,how='outer',on=['filename'])

g1 = pd.merge(f1,f2,how='outer',on=['filename'])
g2 = pd.merge(g1,f3,how='outer',on=['filename'])

cols = ['Mass density, kg/m3', 'Mass density, kg/m3_binary', 'Molar enthalpy of vaporization or sublimation, kJ/mol', 'Excess molar enthalpy (molar enthalpy of mixing), kJ/mol', 'Molar heat capacity at constant pressure, J/K/mol', 'Excess molar heat capacity, J/K/mol']

g2.loc[g2['Mass density, kg/m3'].isnull(), ['Mass density, kg/m3']] = g2.loc[g2['Mass density, kg/m3'].isnull(), 'Mass density, kg/m3'].apply(lambda x: [])
g2.loc[g2['Mass density, kg/m3_binary'].isnull(), ['Mass density, kg/m3_binary']] = g2.loc[g2['Mass density, kg/m3_binary'].isnull(), 'Mass density, kg/m3_binary'].apply(lambda x: [])
g2.loc[g2['Molar enthalpy of vaporization or sublimation, kJ/mol'].isnull(), ['Molar enthalpy of vaporization or sublimation, kJ/mol']] = g2.loc[g2['Molar enthalpy of vaporization or sublimation, kJ/mol'].isnull(), 'Molar enthalpy of vaporization or sublimation, kJ/mol'].apply(lambda x: [])
g2.loc[g2['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'].isnull(), ['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol']] = g2.loc[g2['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'].isnull(), 'Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'].apply(lambda x: [])
g2.loc[g2['Molar heat capacity at constant pressure, J/K/mol'].isnull(), ['Molar heat capacity at constant pressure, J/K/mol']] = g2.loc[g2['Molar heat capacity at constant pressure, J/K/mol'].isnull(), 'Molar heat capacity at constant pressure, J/K/mol'].apply(lambda x: [])
g2.loc[g2['Excess molar heat capacity, J/K/mol'].isnull(), ['Excess molar heat capacity, J/K/mol']] = g2.loc[g2['Excess molar heat capacity, J/K/mol'].isnull(), 'Excess molar heat capacity, J/K/mol'].apply(lambda x: [])

a = g2['Mass density, kg/m3_binary'].tolist() 
b = g2['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'].tolist()
c = g2['Excess molar heat capacity, J/K/mol'].tolist()

a1 = list()
b1 = list()
c1 = list()

for line in a:
    line = list(line)
    line = [i.split('__', 1) for i in line]
    a1.append(line)

for line in b:
    line = list(line)
    line = [i.split('__', 1) for i in line]
    b1.append(line)

for line in c:
    line = list(line)
    line = [i.split('__', 1) for i in line]
    c1.append(line)

g2 = g2.drop(g2.columns[[2,4,6]], axis=1)

g2['Mass density, kg/m3_binary'] = pd.Series(a1, index=g2.index)
g2['Excess molar enthalpy(molar enthalpy of mixing), kJ/mol'] = pd.Series(b1, index=g2.index)
g2['Excess molar heat capacity, J/K/mol'] = pd.Series(c1, index=g2.index)

g2.to_json('journals_for_Ken.json', orient='index')

# This was my attempt for formattting the JSON, but there are certain data types in the dataframe that aren't json serializable and I'm unsure of how to handle that
results = {}
for key, df_gb in g2.groupby('filename'):
    results[str(key)] = df_gb.to_dict('records')


#import json
#print json.dumps(results, indent=4)


