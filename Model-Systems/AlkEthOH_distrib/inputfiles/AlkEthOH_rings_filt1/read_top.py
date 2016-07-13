from parmed.amber import *
import numpy as np
import glob
import pandas as pd

files = glob.glob('./AlkEthOH_r47*.top')

def drop(mylist, m, n):
	mylist = list(mylist)
	del mylist[m::n]	
	return mylist

# Reading in and cleaning up atoms involved in bonds
lst0name = []
lst0 = []
lst00 = []
print("PRINTING BOND PAIRS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	a1 = fin.parm_data['BONDS_INC_HYDROGEN']
	a2 = fin.parm_data['BONDS_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the bond length
	a1 = drop(a1,2,3)
	a2 = drop(a2,2,3)
	
	# Don't need to distinguish between bonded to H or not
	a1.extend(a2)
	
	# Return true atom numbers based on AMBER documentation
	a1 = np.array(a1)/3 + 1
	
	# Subdivide array into those of length 2 to make assigning column titles easier later
	a1 = np.array_split(a1, len(a1)/2)
        a2 = list()
	for line in a1:
		line = list(line)
		a2.append(line)
	
	# Need to create multiple lists for this to work
	# lst0name and lst0 allow me to keep the bond pairs indexed with the molecule
	# lst00 will allow me to create the column names after finding the unique pairs 
	lst0name.append(FileName)
	lst0.append(a1)
	lst00.extend(a2)

# Fun little one liner to get unique pais from lst00
cols0 = [list(x) for x in set(tuple(x) for  x in lst00)]

# Reading in and cleaning up atoms involved in angles
lst1name = []
lst1 = []
lst11 = []
print("PRINTING ANGLE TRIPLETS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	b1 = fin.parm_data['ANGLES_INC_HYDROGEN']
	b2 = fin.parm_data['ANGLES_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the angles
	b1 = drop(b1,3,4)
	b2 = drop(b2,3,4)
	
	# Don't need to distinguish between angles including H or not
	b1.extend(b2)
	
	# Return true atom numbers based on AMBER documentation
	b1 = np.array(b1)/3 + 1
		
	# Subdivide array into those of length 3 to make assigning column titles easier later
	b1 = np.array_split(b1, len(b1)/3)
	b2 = list()
	for line in b1:
		line = list(line)
		b2.append(line)
	
	# Need to create multiple lists for this to work
	# lst1name and lst1 allow me to keep the angle trips indexed with the molecule
	# lst11 will allow me to create the column names after finding the unique trios
	lst1name.append(FileName)
	lst1.append(b1)
	lst11.extend(b2)

# Fun little one liner to get unique trios from lst11
cols1 = [list(x) for x in set(tuple(x) for x in lst11)]

# Reading in and cleaning up atoms involved in dihedrals
lst2 = [] 
lst2name = []
lst22 = []
print("PRINTING DIHEDRAL QUARTETS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	c1 = fin.parm_data['DIHEDRALS_INC_HYDROGEN']
	c2 = fin.parm_data['DIHEDRALS_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the torsions
	c1 = drop(c1,4,5)
	c2 = drop(c2,4,5)
	
	# Don't need to distinguish between torsions including H or not
	c1.extend(c2)
	
	# Return true atom numbers based on AMBER documentation
	for i in range(len(c1)):
		if c1[i] >= 0:
			c1[i] = np.array(c1[i])/3 + 1
		else:
			c1[i] = -(abs(np.array(c1[i]))/3 + 1)
				
	# Subdivide array into those of length 4 to make assigning column titles easier later
	c1 = np.array_split(c1, len(c1)/4)
	c2 = list()
	for line in c1:
		line = list(line)
		c2.append(line)	
	
	# Need to create multiple list for this to work	
	lst2name.append(FileName)
	lst2.append(c1)
	lst22.extend(c2)

# Fun little one liner to get unique trios from lst22
cols2 = [list(x) for x in set(tuple(x) for x in lst22)]


df0 = pd.DataFrame(index = lst0name, columns = cols0)

df0.to_csv("check.csv")
