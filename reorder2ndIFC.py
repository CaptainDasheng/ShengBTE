#!/usr/bin/env python
#

import os
import sys
import numpy as np
import math
import glob


def main(argv = None):

    ifc_file = open('FORCE_CONSTANTS', 'r')
    ifc_data = ifc_file.readlines()
    ifc_file.close()
    ifc_for_sheng=dict()
    
   # f = open('FORCE_CONSTANTS_2ND_reordered', 'w')
   # f_struct = '%6s' %str(ifc_data[0])
    num_pairs =   '%6s' %str(ifc_data[0])
    for index, line in enumerate(ifc_data):
      
      if len(line.split()) < 3 and len(line.split()) > 1:
        dictionary = []
        a = int(line.split()[0])
        b = int(line.split()[1])
        if a > 500 and a < 1001:
          a += 500
	  #print 'a between 500 and 1000', a
          
        elif a > 1000:
          a -= 500
          #print 'a larger than 1000', a
        if b > 500 and b < 1001:
          b += 500
	  #print 'b between 500 and 1000', a
          
        elif b > 1000:
          b -= 500
          #print 'b > 1000:', b
       # f_struct += '%6s' %str(a) + '%6s' %str(b) + '\n' 
        ion_couple = '%6s' %a + '%6s' %b
       # print ion_couple
        starting_line = index + 1 
        for i in range(3):
      	# f_struct += ifc_data[starting_line] 
         dictionary += [[float(ifc_data[starting_line].split()[j]) for j in range (3)]]
         #print dictionary
         #ifc_for_sheng[ion_couple].append([ifc_data[starting_line]])
         starting_line +=1
        ifc_for_sheng[ion_couple] = dictionary
   # print ifc_for_sheng
   # f.write(f_struct)
   # f.close() 

 
    f = open('FORCE_CONSTANTS_reordered', 'w') 
    f_file = '%6s' %num_pairs
    for key1 in sorted(ifc_for_sheng):
       f_file += str(key1) + '\n'
       for m in range(3):
         for l in range (3):
           f_file += '%22s' %str("{0:.15f}".format(ifc_for_sheng[key1][m][l])) 
         f_file += '\n'
    f.write(f_file)
    f.close()   


if __name__ == "__main__":
    import sys
    sys.exit(main())
