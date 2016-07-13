import os
import glob

list_of_files = glob.glob('/inputfiles/*.top')

for fileName in list_of_files:
	fin = open(fileName, 'r')
	data_list = fin.readlines()
	fin.close() # close file

	del data_list[0:17
