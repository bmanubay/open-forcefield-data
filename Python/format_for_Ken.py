import thermopyl as th 
from thermopyl import thermoml_lib
import cirpy
import numpy as np
import pandas as pd
from sklearn.externals.joblib import Memory

#def format_for_Ken(filenames):
    # filenames and properties are both just arrays of strings to use as filter keys for the function
    # import ThermoML as pandas df
#    df = th.pandas_dataframe()
    
    # strip leading and trailing strings off of filename info
#    df['filename'] = df['filename'].map(lambda x: x.lstrip('/home/bmanubay/.thermoml/'))
#    df['filename'] = df['filename'].map(lambda x: x.replace(' ', '')[:-4])
    
    # keep only data points in filenames parameter
#    df = df[df.filename.isin(filenames)]

    # list unique components for each journal article in filenames
#    df1 = df.groupby('filename').components.unique() 

    # get data point counts for each of the properties I'm interested in and drop columns with NaN or 0
    # merge with df1 on filename
    # store rows as array
    # save files as json

#    return df1

#filenames = ["je900468u"]

#a = format_for_Ken(filenames)

import numpy as np
import pandas as pd

# Read in property data produced from puresrch.py and binsrch.py
pathdfp = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Property data/"
pathdfm = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Binary-Mixtures/Property data/"

a1 = pd.read_csv(pathdfp+"dens_pure.csv", sep=';')   
a2 = pd.read_csv(pathdfp+"hvap_pure.csv", sep=';')
a3 = pd.read_csv(pathdfp+"cpmol_pure.csv", sep=';')

b1 = pd.read_csv(pathdfm+"dens_bin.csv", sep=';')
b2 = pd.read_csv(pathdfm+"eme_bin.csv", sep=';')
b3 = pd.read_csv(pathdfm+"emcp_bin.csv", sep=';')

# Create new columns for combined SMILES strings on mixtures
b1 = b1.rename(columns = {'components':'components_mix'})
b2 = b2.rename(columns = {'components':'components_mix'})
b3 = b3.rename(columns = {'components':'components_mix'})

# Merge property sets using filename as index
c1 = pd.merge(a1,b1,how='outer',on=['filename'])
c2 = pd.merge(a2,b2,how='outer',on=['filename'])
c3 = pd.merge(a3,b3,how='outer',on=['filename'])

# Find unique compoounds/mixtures for each journal in merged datasets 
d1 = c1.groupby('filename').components.unique()
d1 = d1.reset_index()
d1.columns = ['filename', 'Mass density, kg/m3']
e1 = c1.groupby('filename').components_mix.unique()
e1 = e1.reset_index()
e1.columns = ['filename', 'Mass density, kg/m3_binary']
d2 = c2.groupby('filename').components.unique()
d2 = d2.reset_index()
d2.columns = ['filename', 'Molar enthalpy of vaporization or sublimation, kJ/mol']
e2 = c2.groupby('filename').components_mix.unique()
e2 = e2.reset_index()
e2.columns = ['filename', 'Excess molar enthalpy (molar enthalpy of mixing), kJ/mol']
d3 = c3.groupby('filename').components.unique()
d3 = d3.reset_index()
d3.columns = ['filename', 'Molar heat capacity at constant pressure, J/K/mol']
e3 = c3.groupby('filename').components_mix.unique()
e3 = e3.reset_index()
e3.columns = ['filename', 'Excess molar heat capacity, J/K/mol']

# Merge the unique component/mixture list with filename as index
f1 = pd.merge(d1,e1,how='outer',on=['filename'])
f2 = pd.merge(d2,e2,how='outer',on=['filename'])
f3 = pd.merge(d3,e3,how='outer',on=['filename'])

g1 = pd.merge(f1,f2,how='outer',on=['filename'])
g2 = pd.merge(g1,f3,how='outer',on=['filename'])

# Convert data type to string to allow simpler application of counting process
g2 = g2.astype('str')

# Format strings in columns
g2['Mass density, kg/m3'] = g2['Mass density, kg/m3'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Mass density, kg/m3_binary'] = g2['Mass density, kg/m3_binary'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Molar enthalpy of vaporization or submlimation, kJ/mol'] = g2['Molar enthalpy of vaporization or sublimation, kJ/mol'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'] = g2['Excess molar enthalpy (molar enthalpy of mixing), kJ/mol'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Molar heat capacity at constant pressure, J/K/mol'] = g2['Molar heat capacity at constant pressure, J/K/mol'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Excess molar heat capacity, J/K/mol'] = g2['Excess molar heat capacity, J/K/mol'].map(lambda x: x.lstrip('"[').rstrip(']"'))


g2.to_json('journals_for_Ken.json', orient='index')
