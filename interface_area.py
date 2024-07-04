## run code within pymol.. while you have all files and codes in the same directory.
## i suggest to first open any PDB file in the same folder to set the default path to current directory. and then run the pyhton filename in the pymol command line. It will run for all PDBs in the directory and give you CSV file

import os
import pymol
dir =  os.getcwd()
files = [x for x in os.listdir(dir) if x.endswith(".pdb")]
file3 = open("int_area.csv","w")

for i in files:
	cmd.set("dot_solvent", "on")
	cmd.set("dot_density",4)
	cmd.set("solvent_radius", 1.4)
	cmd.delete("all")
	cmd.load(i)
	ALL = cmd.get_area("all")
	cmd.remove("chain A")
	HL = cmd.get_area("chain R")
	cmd.delete("all")
	cmd.load(i)
	cmd.remove("chain R")
	E = cmd.get_area("chain A")
	cmd.delete("all")
	
	val = HL+E-ALL
	file3.write(i+","+str(HL)+","+str(E)+","+str(ALL)+","+str(val/2)+"\n")
	print (i)
	print (HL)
	print (E)
	print (ALL)
file3.close()
