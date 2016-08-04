# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 17:38:07 2016

@author: Bryce Manubay
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

purepath = "C:\\Users\\Bryce Manubay\\Desktop\\UCB\\Shirts Group\\open-forcefield-data-master\\Pure-Solvents\\Component counts\\"
mixdpath = "C:\\Users\\Bryce Manubay\\Desktop\\UCB\\Shirts Group\\open-forcefield-data-master\\Binary-Mixtures\\Component counts\\"

a0 = pd.read_csv(purepath+"purecomp_counts_all.csv", sep=";", usecols=["Count", "SMILES"])
a0 = a0.rename(columns={"Count":"PureAll"})
a1 = pd.read_csv(purepath+"purecomp_counts_dens.csv", sep=";", usecols=["Count", "SMILES"])
a1 = a1.rename(columns={"Count":"PureDens"})
a2 = pd.read_csv(purepath+"purecomp_counts_sos.csv", sep=";", usecols=["Count", "SMILES"])
a2 = a2.rename(columns={"Count":"PureSOS"})
a3 = pd.read_csv(purepath+"purecomp_counts_dielec.csv", sep=";", usecols=["Count", "SMILES"])
a3 = a3.rename(columns={"Count":"PureDielec"})
a4 = pd.read_csv(purepath+"purecomp_counts_cpmol.csv", sep=";", usecols=["Count", "SMILES"])
a4 = a4.rename(columns={"Count":"PureCP"})
a5 = pd.read_csv(purepath+"purecomp_counts_hvap.csv", sep=";", usecols=["Count", "SMILES"])
a5 = a5.rename(columns={"Count":"PureHvap"})

aa0 = pd.read_csv(purepath+"purecomp_counts_all_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa0 = aa0.rename(columns={"Count":"PureAllAlk"})
aa1 = pd.read_csv(purepath+"purecomp_counts_dens_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa1 = aa1.rename(columns={"Count":"PureDensAlk"})
aa2 = pd.read_csv(purepath+"purecomp_counts_sos_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa2 = aa2.rename(columns={"Count":"PureSOSAlk"})
aa3 = pd.read_csv(purepath+"purecomp_counts_dielec_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa3 = aa3.rename(columns={"Count":"PureDielecAlk"})
aa4 = pd.read_csv(purepath+"purecomp_counts_cpmol_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa4 = aa4.rename(columns={"Count":"PureCPAlk"})
aa5 = pd.read_csv(purepath+"purecomp_counts_hvap_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa5 = aa5.rename(columns={"Count":"PureHvapAlk"})

b0 = pd.read_csv(mixdpath+"bincomp_counts_all.csv", sep=";", usecols=["Count", "SMILES"])
b0 = b0.rename(columns={"Count":"BinAll"})
b1 = pd.read_csv(mixdpath+"bincomp_counts_dens.csv", sep=";", usecols=["Count", "SMILES"])
b1 = b1.rename(columns={"Count":"BinDens"})
b2 = pd.read_csv(mixdpath+"bincomp_counts_sos.csv", sep=";", usecols=["Count", "SMILES"])
b2 = b2.rename(columns={"Count":"BinSOS"})
b3 = pd.read_csv(mixdpath+"bincomp_counts_dielec.csv", sep=";", usecols=["Count", "SMILES"])
b3 = b3.rename(columns={"Count":"BinDielec"})
b4 = pd.read_csv(mixdpath+"bincomp_counts_eme.csv", sep=";", usecols=["Count", "SMILES"])
b4 = b4.rename(columns={"Count":"BinEME"})
b5 = pd.read_csv(mixdpath+"bincomp_counts_emcp.csv", sep=";", usecols=["Count", "SMILES"])
b5 = b5.rename(columns={"Count":"BinEMC"})
b6 = pd.read_csv(mixdpath+"bincomp_counts_emv.csv", sep=";", usecols=["Count", "SMILES"])
b6 = b6.rename(columns={"Count":"BinEMV"})
b7 = pd.read_csv(mixdpath+"bincomp_counts_actcoeff.csv", sep=";", usecols=["Count", "SMILES"])
b7 = b7.rename(columns={"Count":"BinAct"})
bb0 = pd.read_csv(mixdpath+"bincomp_counts_all_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb0 = bb0.rename(columns={"Count":"BinAllAlk"})
bb1 = pd.read_csv(mixdpath+"bincomp_counts_dens_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb1 = bb1.rename(columns={"Count":"BinDensAlk"})
bb2 = pd.read_csv(mixdpath+"bincomp_counts_sos_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb2 = bb2.rename(columns={"Count":"BinSOSAlk"})
bb3 = pd.read_csv(mixdpath+"bincomp_counts_dielec_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb3 = bb3.rename(columns={"Count":"BinDielecAlk"})
bb4 = pd.read_csv(mixdpath+"bincomp_counts_eme_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb4 = bb4.rename(columns={"Count":"BinEMEAlk"})
bb5 = pd.read_csv(mixdpath+"bincomp_counts_emcp_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb5 = bb5.rename(columns={"Count":"BinEMCAlk"})
bb6 = pd.read_csv(mixdpath+"bincomp_counts_emv_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb6 = bb6.rename(columns={"Count":"BinEMVAlk"})
bb7 = pd.read_csv(mixdpath+"bincomp_counts_actcoeff_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb7 = bb7.rename(columns={"Count":"BinActAlk"})

# Merge dataframes on SMILES column
dfpur = a0.merge(a1,how='outer',on='SMILES').merge(a2,how='outer',on='SMILES').merge(a3,how='outer',on='SMILES').merge(a4,how='outer',on='SMILES').merge(a5,how='outer',on='SMILES')
dfpur['PureScore'] = dfpur.notnull().sum(axis=1) - 2
dfbin = b0.merge(b1,how='outer',on='SMILES').merge(b2,how='outer',on='SMILES').merge(b3,how='outer',on='SMILES').merge(b4,how='outer',on='SMILES').merge(b5,how='outer',on='SMILES').merge(b6,how='outer',on='SMILES').merge(b7,how='outer',on='SMILES')
dfbin['BinScore'] = dfbin.notnull().sum(axis=1) - 2
df = dfpur.merge(dfbin,how='outer',on='SMILES')

dfpurAlk = aa0.merge(aa1,how='outer',on='SMILES').merge(aa2,how='outer',on='SMILES').merge(aa3,how='outer',on='SMILES').merge(aa4,how='outer',on='SMILES').merge(aa5,how='outer',on='SMILES')
dfpurAlk['PureScore'] = dfpurAlk.notnull().sum(axis=1) - 2
dfbinAlk = bb0.merge(bb1,how='outer',on='SMILES').merge(bb2,how='outer',on='SMILES').merge(bb3,how='outer',on='SMILES').merge(bb4,how='outer',on='SMILES').merge(bb5,how='outer',on='SMILES').merge(bb6,how='outer',on='SMILES').merge(bb7,how='outer',on='SMILES')
dfbinAlk['BinScore'] = dfbinAlk.notnull().sum(axis=1) - 2
dff = dfpurAlk.merge(dfbinAlk,how='outer',on='SMILES')

ind1 = (df.PureScore>=1) & (df.BinScore>=1)
df1 = df[ind1] 
len11 = len(df1.SMILES)
ind2 = (df.PureScore>=2) & (df.BinScore>=2)
df2 = df[ind2] 
len22 = len(df2.SMILES)
ind3 = (df.PureScore>=3) & (df.BinScore>=3)
df3 = df[ind3] 
len33 = len(df3.SMILES)
ind4 = (df.PureScore>=4) & (df.BinScore>=4)
df4 = df[ind4] 
len44 = len(df4.SMILES)
ind5 = (df.PureScore>=5) & (df.BinScore>=5)
df5 = df[ind5] 
len55 = len(df5.SMILES)
ind6 = (df.PureScore>=5) & (df.BinScore>=6)
df6 = df[ind6] 
len56 = len(df6.SMILES)
ind7 = (df.PureScore>=5) & (df.BinScore>=7)
df7 = df[ind7] 
len57 = len(df7.SMILES)

ind11 = (dff.PureScore>=1) & (dff.BinScore>=1)
df11 = dff[ind11] 
len11Alk = len(df11.SMILES)
ind22 = (dff.PureScore>=2) & (dff.BinScore>=2)
df22 = dff[ind22] 
len22Alk = len(df22.SMILES)
ind33 = (dff.PureScore>=3) & (dff.BinScore>=3)
df33 = dff[ind33] 
len33Alk = len(df33.SMILES)
ind44 = (dff.PureScore>=4) & (dff.BinScore>=4)
df44 = dff[ind44] 
len44Alk = len(df44.SMILES)
ind55 = (dff.PureScore>=5) & (dff.BinScore>=5)
df55 = dff[ind55] 
len55Alk = len(df55.SMILES)
ind66 = (dff.PureScore>=5) & (dff.BinScore>=6)
df66 = dff[ind66] 
len56Alk = len(df66.SMILES)
ind77 = (dff.PureScore>=5) & (dff.BinScore>=7)
df77 = dff[ind77] 
len57Alk = len(df77.SMILES)

# Number of molecules in each property set
purecountall = len(a0.SMILES)
purecountden = len(a1.SMILES)
purecountsos = len(a2.SMILES)
purecountdie = len(a3.SMILES)
purecountcpm = len(a4.SMILES)
purecounthva = len(a5.SMILES)
purecountallAlk = len(aa0.SMILES)
purecountdenAlk = len(aa1.SMILES)
purecountsosAlk = len(aa2.SMILES)
purecountdieAlk = len(aa3.SMILES)
purecountcpmAlk = len(aa4.SMILES)
purecounthvaAlk = len(aa5.SMILES)
bincountall = len(b0.SMILES)
bincountden = len(b1.SMILES)
bincountsos = len(b2.SMILES)
bincountdie = len(b3.SMILES)
bincounteme = len(b4.SMILES)
bincountemc = len(b5.SMILES)
bincountemv = len(b6.SMILES)
bincountact = len(b7.SMILES)
bincountallAlk = len(bb0.SMILES)
bincountdenAlk = len(bb1.SMILES)
bincountsosAlk = len(bb2.SMILES)
bincountdieAlk = len(bb3.SMILES)
bincountemeAlk = len(bb4.SMILES)
bincountemcAlk = len(bb5.SMILES)
bincountemvAlk = len(bb6.SMILES)
bincountactAlk = len(bb7.SMILES)

def autolabel(rects, ax):
    # Get y-axis height to calculate label position from.
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom

    for rect in rects:
        height = rect.get_height()

        # Fraction of axis height taken up by this rectangle
        p_height = (height / y_height)

        # If we can fit the label above the column, do that;
        # otherwise, put it inside the column.
        if p_height > 0.95: # arbitrary; 95% looked good to me.
            label_position = height - (y_height * 0.05)
        else:
            label_position = height + (y_height * 0.01)

        ax.text(rect.get_x() + rect.get_width()/2., label_position,
                '%d' % int(height),
                ha='center', va='bottom')

# make plots
MolDat = [purecountall, purecountden, purecountsos, purecountdie, purecountcpm, purecounthva, bincountall, bincountden, bincountsos, bincountdie, bincounteme, bincountemc, bincountemv, bincountact]
MolDatAlk = [purecountallAlk, purecountdenAlk, purecountsosAlk, purecountdieAlk, purecountcpmAlk, purecounthvaAlk, bincountallAlk, bincountdenAlk, bincountsosAlk, bincountdieAlk, bincountemeAlk, bincountemcAlk, bincountemvAlk, bincountactAlk]

N = 14
ind = np.arange(N)  # the x locations for the groups
width = 0.2      # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, MolDat, width, color='r')

rects2 = ax.bar(ind + width, MolDatAlk, width, color='y')

autolabel(rects1, ax)
autolabel(rects2, ax)

ax.set_ylabel('Molecule Count')
ax.set_title('Unique Molecules per Property of Interest')
ax.set_xticks(ind + width)
ax.set_xticklabels(('Pure All', 'Pure Density', 'Pure Speed of Sound', 'Pure Dielectric Constant', 'Pure Heat Capacity', 'Pure Enthalpy of Vaporization', 'Binary All', 'Binary Density', 'Binary Speed of Sound', 'Binary Dielectric Constant', 'Binary Excess Molar Enthalpy', 'Binary Excess Molar Heat Capacity', 'Binary Excess Molar Volume', 'Binary Activity Coefficient'), rotation=90)

ax.legend((rects1[0], rects2[0]), ('All data', 'AlkEthOH set'))

plt.show()

SpreadDat = [len11,len22,len33,len44,len55,len56,len57]
SpreadDatAlk = [len11Alk,len22Alk,len33Alk,len44Alk,len55Alk,len56Alk,len57Alk]

N = 7
ind = np.arange(N)  # the x locations for the groups
width = 0.2      # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, SpreadDat, width, color='r')

rects2 = ax.bar(ind + width, SpreadDatAlk, width, color='y')

autolabel(rects1, ax)
autolabel(rects2, ax)

ax.set_ylabel('Molecule Count')
ax.set_title('Coupled Pure and Binary Property Coverage')
ax.set_xticks(ind + width)
ax.set_xticklabels(('>=1 Pure and >=1 Binary', '>=2 Pure and >=2 Binary', '>=3 Pure and >=3 Binary', '>=4 Pure and >=4 Binary', '>=5 Pure and >=5 Binary', '>=5 Pure and >=6 Binary', '>=5 Pure and >=7 Binary'), rotation=90)

ax.legend((rects1[0], rects2[0]), ('All data', 'AlkEthOH set'))

plt.show()
