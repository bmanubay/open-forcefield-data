import thermopyl as th
import pandas as pd
from pdb import set_trace


df = th.pandas_dataframe()

experiments = ['Molar volume, m3/mol', 'Molar volume, m3/mol_std', 'Molar heat capacity at constant pressure, J/K/mol', 'Molar heat capacity at constant pressure, J/K/mol_std']

ind_list = [df[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
df = df.ix[ind]

columns_keep = ['Molar volume, m3/mol', 'Molar volume, m3/mol_std', 'Molar heat capacity at constant pressure, J/K/mol', 'Molar heat capacity at constant pressure, J/K/mol_std', 'Temperature, K', 'Pressure, kPa', 'components', 'filename', 'phase']

df = df[df['components']=='cyclohexane']
df = df[df['Temperature, K']==293.15]
df = df[df['Pressure, kPa']==101.325]

df = df[columns_keep]

print df

df.to_csv('tabular_cychex_data_for_Josh.csv',sep=';')

