import thermopyl as th
import pandas as pd
from pdb import set_trace


df = th.pandas_dataframe()
datatypes = list(df.columns)
print datatypes

file_name = open('DOIs-1.txt', 'r')
DOI = file_name.readlines()


# Remove first part of DOI before '/'
DOI_trunc = []
for i in DOI:
    try:
        i = i.split('/')[1]
        i = i.rstrip('\n')
        DOI_trunc.append(i)
    except IndexError:
        print "Invalid file name" 
        continue 

# Prepend path and append '.xml' to each entry in DOI
path = '/home/brma3379/.thermoml/'
DOI = [path+i+'.xml' for i in DOI_trunc]

# Slice data frame on DOI list
df = df[df.filename.isin(DOI)]

experiments = ["Mass density, kg/m3","Molar enthalpy of vaporization or sublimation, kJ/mol","Excess molar enthalpy (molar enthalpy of mixing), kJ/mol"]

ind_list = [df[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
df = df.ix[ind]

columns_keep = ["Mass density, kg/m3","Molar enthalpy of vaporization or sublimation, kJ/mol","Excess molar enthalpy (molar enthalpy of mixing), kJ/mol","Temperature, K","Pressure, kPa","filename","components","phase","Mole fraction"]

df = df[columns_keep]

# Extract rows with temperature between 128 and 399 K
df = df[df['Temperature, K'] > 273.]
df = df[df['Temperature, K'] < 373.]

# Extract rows with pressure between 101.325 kPa and 101325 kPa
df = df[df['Pressure, kPa'] > 50.]
df = df[df['Pressure, kPa'] < 203.]

print df

df.to_csv('tabular_data_for_MKG.csv',sep=';')

