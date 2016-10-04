# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 17:38:07 2016

@author: Bryce Manubay
"""
import matplotlib as mpl

mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb
import sys

purepath = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Pure-Solvents/Component counts/"
mixdpath = "/home/bmanubay/.thermoml/tables/Ken/open-forcefield-data/Binary-Mixtures/Component counts/"

a0 = pd.read_csv(purepath+"purecomp_counts_all_full.csv", sep=";", usecols=["Count", "SMILES"])
a0 = a0.rename(columns={"Count":"PureAll"})
a1 = pd.read_csv(purepath+"purecomp_counts_dens_full.csv", sep=";", usecols=["Count", "SMILES"])
a1 = a1.rename(columns={"Count":"PureDens"})
a2 = pd.read_csv(purepath+"purecomp_counts_sos_full.csv", sep=";", usecols=["Count", "SMILES"])
a2 = a2.rename(columns={"Count":"PureSOS"})
a3 = pd.read_csv(purepath+"purecomp_counts_dielec_full.csv", sep=";", usecols=["Count", "SMILES"])
a3 = a3.rename(columns={"Count":"PureDielec"})
a4 = pd.read_csv(purepath+"purecomp_counts_cpmol_full.csv", sep=";", usecols=["Count", "SMILES"])
a4 = a4.rename(columns={"Count":"PureCP"})
a5 = pd.read_csv(purepath+"purecomp_counts_hvap_full.csv", sep=";", usecols=["Count", "SMILES"])
a5 = a5.rename(columns={"Count":"PureHvap"})

aa0 = pd.read_csv(purepath+"purecomp_counts_all_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa0 = aa0.rename(columns={"Count":"PureAllAlk"})
aa1 = pd.read_csv(purepath+"purecomp_counts_dens_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa1 = aa1.rename(columns={"Count":"PureDensAlk"})
aa2 = pd.read_csv(purepath+"purecomp_counts_sos_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa2 = aa2.rename(columns={"Count":"PureSOSAlk"})
aa3 = pd.read_csv(purepath+"purecomp_counts_dielec_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa3 = aa3.rename(columns={"Count":"PureDielecAlk"})
aa4 = pd.read_csv(purepath+"purecomp_counts_cpmol_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa4 = aa4.rename(columns={"Count":"PureCPAlk"})
aa5 = pd.read_csv(purepath+"purecomp_counts_hvap_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
aa5 = aa5.rename(columns={"Count":"PureHvapAlk"})

b0 = pd.read_csv(mixdpath+"bincomp_counts_all_full.csv", sep=";", usecols=["Count", "SMILES"])
b0 = b0.rename(columns={"Count":"BinAll"})
b1 = pd.read_csv(mixdpath+"bincomp_counts_dens_full.csv", sep=";", usecols=["Count", "SMILES"])
b1 = b1.rename(columns={"Count":"BinDens"})
b2 = pd.read_csv(mixdpath+"bincomp_counts_sos_full.csv", sep=";", usecols=["Count", "SMILES"])
b2 = b2.rename(columns={"Count":"BinSOS"})
b3 = pd.read_csv(mixdpath+"bincomp_counts_dielec_full.csv", sep=";", usecols=["Count", "SMILES"])
b3 = b3.rename(columns={"Count":"BinDielec"})
b4 = pd.read_csv(mixdpath+"bincomp_counts_eme_full.csv", sep=";", usecols=["Count", "SMILES"])
b4 = b4.rename(columns={"Count":"BinEME"})
b5 = pd.read_csv(mixdpath+"bincomp_counts_emcp_full.csv", sep=";", usecols=["Count", "SMILES"])
b5 = b5.rename(columns={"Count":"BinEMC"})
b6 = pd.read_csv(mixdpath+"bincomp_counts_emv_full.csv", sep=";", usecols=["Count", "SMILES"])
b6 = b6.rename(columns={"Count":"BinEMV"})
b7 = pd.read_csv(mixdpath+"bincomp_counts_actcoeff_full.csv", sep=";", usecols=["Count", "SMILES"])
b7 = b7.rename(columns={"Count":"BinAct"})
bb0 = pd.read_csv(mixdpath+"bincomp_counts_all_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb0 = bb0.rename(columns={"Count":"BinAllAlk"})
bb1 = pd.read_csv(mixdpath+"bincomp_counts_dens_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb1 = bb1.rename(columns={"Count":"BinDensAlk"})
bb2 = pd.read_csv(mixdpath+"bincomp_counts_sos_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb2 = bb2.rename(columns={"Count":"BinSOSAlk"})
bb3 = pd.read_csv(mixdpath+"bincomp_counts_dielec_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb3 = bb3.rename(columns={"Count":"BinDielecAlk"})
bb4 = pd.read_csv(mixdpath+"bincomp_counts_eme_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb4 = bb4.rename(columns={"Count":"BinEMEAlk"})
bb5 = pd.read_csv(mixdpath+"bincomp_counts_emcp_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb5 = bb5.rename(columns={"Count":"BinEMCAlk"})
bb6 = pd.read_csv(mixdpath+"bincomp_counts_emv_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
bb6 = bb6.rename(columns={"Count":"BinEMVAlk"})
bb7 = pd.read_csv(mixdpath+"bincomp_counts_actcoeff_full_AlkEthOH.csv", sep=";", usecols=["Count", "SMILES"])
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

# measure atom type diversity in AlkEthOH filter set 
aa0alkane = aa0[aa0.SMILES.str.contains('O') == False]
bb0alkane = bb0[bb0.SMILES.str.contains('O') == False]

aa0alco = aa0[aa0.SMILES.str.contains('CO') == True]
aa0alco['if_o'] = aa0.SMILES.str[-1:]
aa0alco = aa0alco[aa0alco.if_o.str.contains('O') == True]
aa0alco = aa0alco.drop('if_o', 1)
bb0alco = bb0[bb0.SMILES.str.contains('CO') == True]
bb0alco['if_o'] = bb0.SMILES.str[-1:]
bb0alco = bb0alco[bb0alco.if_o.str.contains('O') == True]
bb0alco = bb0alco.drop('if_o', 1)

aa0ether = aa0[(aa0.SMILES.str.contains('COC') == True) | (aa0.SMILES.str.contains('(C)OC') == True) | (aa0.SMILES.str.contains('CO(C)') == True)] 
bb0ether = bb0[(bb0.SMILES.str.contains('COC') == True) | (bb0.SMILES.str.contains('(C)OC') == True) | (bb0.SMILES.str.contains('CO(C)') == True)] 

puralkpts = aa0alkane.PureAllAlk.sum()
puralcpts = aa0alco.PureAllAlk.sum()
purethpts = aa0ether.PureAllAlk.sum()

binalkpts = bb0alkane.BinAllAlk.sum()
binalcpts = bb0alco.BinAllAlk.sum()
binethpts = bb0ether.BinAllAlk.sum()

puralktyp = aa0alkane.PureAllAlk.count()
puralctyp = aa0alco.PureAllAlk.count()
purethtyp = aa0ether.PureAllAlk.count()

binalktyp = bb0alkane.BinAllAlk.count()
binalctyp = bb0alco.BinAllAlk.count()
binethtyp = bb0ether.BinAllAlk.count()


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
rects1 = ax.bar(ind, MolDat, width, color='r', align='center')

rects2 = ax.bar(ind + width, MolDatAlk, width, color='y', align='center')

autolabel(rects1, ax)
autolabel(rects2, ax)

ax.set_ylabel('Molecule Count')
ttl = ax.set_title('Unique Molecules per Property of Interest')
ttl.set_position([.5, 1.05])
ax.set_xticks(ind + width)
xlabel = ax.set_xticklabels(('Pure All', 'Pure Density', 'Pure Speed of Sound', 'Pure Dielectric Constant', 'Pure Heat Capacity', 'Pure Enthalpy of Vaporization', 'Binary All', 'Binary Density', 'Binary Speed of Sound', 'Binary Dielectric Constant', 'Binary Excess Molar Enthalpy', 'Binary Excess Molar Heat Capacity', 'Binary Excess Molar Volume', 'Binary Activity Coefficient'), rotation=90)

ax.legend((rects1[0], rects2[0]), ('All data', 'AlkEthOH set'))

plt.savefig('Molecules_per_property.png', bbox_inches='tight')

SpreadDat = [len11,len22,len33,len44,len55,len56,len57]
SpreadDatAlk = [len11Alk,len22Alk,len33Alk,len44Alk,len55Alk,len56Alk,len57Alk]

N = 7
ind = np.arange(N)  # the x locations for the groups
width = 0.2      # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, SpreadDat, width, color='r', align='center')

rects2 = ax.bar(ind + width, SpreadDatAlk, width, color='y', align='center')

autolabel(rects1, ax)
autolabel(rects2, ax)

ax.set_ylabel('Molecule Count')
ttl = ax.set_title('Coupled Pure and Binary Property Coverage')
ttl.set_position([.5, 1.05])
ax.set_xticks(ind + width)
xlabel = ax.set_xticklabels(('>=1 Pure and >=1 Binary', '>=2 Pure and >=2 Binary', '>=3 Pure and >=3 Binary', '>=4 Pure and >=4 Binary', '>=5 Pure and >=5 Binary', '>=5 Pure and >=6 Binary', '>=5 Pure and >=7 Binary'), rotation=90)

ax.legend((rects1[0], rects2[0]), ('All data', 'AlkEthOH set'))

plt.savefig('Coupled_property_coverage.png', bbox_inches='tight')

CountDat = [puralkpts,puralcpts,purethpts,binalkpts,binalcpts,binethpts]
TypeDat = [puralktyp,puralctyp,purethtyp,binalktyp,binalctyp,binethtyp]

N = 6
ind = np.arange(N)
width = 0.2

fig, ax = plt.subplots()

rects1 = ax.bar(ind,CountDat,width,color='r',align='center')

autolabel(rects1,ax)

ax.set_ylabel('Data Point Count')
ttl = ax.set_title('Data point count per chemical environment in AlkEthOH filtered data')
ttl.set_position([.5, 1.05])
ax.set_xticks(ind+width)
xlabel = ax.set_xticklabels(('Alkanes Pure','Alcohols Pure','Ethers Pure','Alkanes Binary','Alcohols Binary','Ethers Binary'), rotation=90)

plt.savefig('data_per_chemistry.png', bbox_inches='tight')

N = 6
ind = np.arange(N)
width = 0.2

fig, ax = plt.subplots()

rects1 = ax.bar(ind,TypeDat,width,color='r',align='center')

autolabel(rects1,ax)

ax.set_ylabel('Molecule Count')
ttl = ax.set_title('Molecules per chemical environment in AlkEthOH filtered data')
ttl.set_position([.5, 1.05])
ax.set_xticks(ind+width)
xlabel = ax.set_xticklabels(('Alkanes Pure','Alcohols Pure','Ethers Pure','Alkanes Binary','Alcohols Binary','Ethers Binary'), rotation=90)

plt.savefig('types_per_chemistry.png', bbox_inches='tight')

