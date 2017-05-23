#!/usr/bin/env python
#

import os
import sys
import numpy as np
import math
import glob


def main(argv = None):

    ifc_file = open('FORCE_CONSTANTS_3RD', 'r')
    ifc_data = ifc_file.readlines()
    ifc_file.close()
    ifc_for_sheng=dict()
 
    num_pairs =   ifc_data[0].split()[0]
    ifc_data[0] +=str(1)
    print num_pairs
    print ifc_data[0]
    for index, line in enumerate(ifc_data):
      
      if len(line.split()) == 1:
	print int(line)        
	cells = []
        starting_line = index + 1 
        for i in range(2):
         cells += [[float(ifc_data[starting_line].split()[j]) for j in range (3)]]
         starting_line +=1
        for k in range(1):
         cells += [[int(ifc_data[starting_line].split()[j]) for j in range (3)]]
         starting_line +=1
        for m in range(27):
         cells += [[float(ifc_data[starting_line].split()[j]) for j in range (4)]]
         starting_line +=1
        ifc_for_sheng[int(line)] = cells
#    print ifc_for_sheng
    # This part changes indices
    for key in ifc_for_sheng:
	for i in range(3):
		atom_index = int(ifc_for_sheng[key][2][i])	
		if atom_index > 4 and atom_index < 9:       
                        atom_index += 4
                        ifc_for_sheng[key][2][i] = atom_index
		elif atom_index > 8:
                        atom_index -= 4
                        ifc_for_sheng[key][2][i] = atom_index 
#    print ifc_for_sheng

    f = open('FORCE_CONSTANTS_3RD_reordered', 'w') 
    f_file = '%-1s' %num_pairs + '\n' + '\n'
    for key1 in sorted(ifc_for_sheng):
       f_file += '%5s' %str(key1) + '\n'
       # cell format
       for m in range(2):
         for l in range (3):
           f_file += '%-18s' %str("{0:.10e}".format(ifc_for_sheng[key1][m][l])) 
         f_file += '\n'
       # atom indices
       for m in range(2,3):
         for l in range (3):
          if l == 0:           
		f_file += '%6s' %str("{0:.0f}".format(ifc_for_sheng[key1][m][l])) 
	  else:
		f_file += '%7s' %str("{0:.0f}".format(ifc_for_sheng[key1][m][l]))
         f_file += '\n'
       # forces format
       for m in range(3,30):
         for l in range (4):
           if l == 3:
           	f_file += '%21s' %str("{0:.10e}".format(ifc_for_sheng[key1][m][l]))
           elif l == 0:
                f_file += '%2s' %str("{0:.0f}".format(ifc_for_sheng[key1][m][l]))
           else:
                f_file += '%3s' %str("{0:.0f}".format(ifc_for_sheng[key1][m][l])) 
         f_file += '\n'
       f_file += '\n'
    f.write(f_file)
    f.close()   


if __name__ == "__main__":
    import sys
    sys.exit(main())
