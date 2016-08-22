#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage ./castep2shengbte.py <seedname>
# where seedname is the name of the .castep file.
#
#
# v0.2 11 May 16 - now loads castep IFC data into a dictionary
# for each supercell. Each supercell contains ions^2 IFC matrices.
#   
# TODO write all permutations between atoms to a file


import os.path
import glob
import os
import sys
import math
import numpy as np

# Convert this to castep. Will use a castep file for the moment

########################################################################
# CASTEP


def main(argv = None):
    print '========================================================='
    print '||              CASTEP 2 ShengBTE Interface            ||'
    print '||                     version 0.2                     ||'
    print '||                     11 May 2016                     ||'
    print '||-----------------------------------------------------||'
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
       # Avoid ugly errors
       print '|| Usage: castep2boltz <seedname> <optional arguments> ||'
       print '|| optional arguments: "so" (for SOC runs) ...         ||'
       print '||          and "down" (for spin down calculations)    ||'
       print '||-----------------------------------------------------||'
       print '|| UNSUCCESSFUL! READ Usage above                      ||'
       print '...'
       sys.exit()
       
    # Define <seedname>
    prefix = argv[1]

    # Help menu, it shows the message and stops the process
    help = ['h', '-h','--h', 'help', '-help', '--help']
    for i in help:
      if i in argv:
            print '|| Usage: castep2boltz <seedname> <optional arguments> ||'
            print '|| optional arguments: "=fgnt6so" (for SOC runs) ...         ||'
            print '||          and "down" (for spin down calculations)    ||'
            print '========================================================='
            sys.exit()


    """
    Return all the relevant information contained in a .castep file.
    """
    castep_file = prefix + '.castep'
    castep_file = open(castep_file, 'r')
    castep_data = castep_file.readlines()
    castep_file.close()

       
    # hardcoded variables which would be extracted from .castep at a later stage.   
    for line in castep_data:
      if 'Total number of ions in cell' in line:
         ions = int(line.split()[7])
      elif 'MP grid size for SCF calculation is' in line:
         print 'Enable later'   
         #nx = int(line.split()[7])
         #ny = int(line.split()[8])
         #nz = int(line.split()[9])
         #supercell_size = nx*ny*nz
    # supercell
    nx=5
    ny=5
    nz=5
    supercell_size = nx*ny*nz
    
    ifc_data=dict()   
    for i in range(supercell_size):
        supercell = "Supercell" + "%4s" %(i+1)
        ifc_data[supercell]=[]         
        for n in range(ions*ions):
           ifc_data[supercell].append(np.empty((3,3)))

        
    # Define a number which gives the number of lines we have per atom per space index in force constants matrix
    l_castep = int(math.ceil(float(ions)/float(2)))
    for index, line in enumerate(castep_data):
         # Extract information for all supercells
         for z in range(supercell_size):
           # supercell_3nx3n_matrix used to avoid castep format for systems with more than 2 atoms in the unit cell
           supercell_3nx3n_matrix = []
           supercell = "Supercell" + "%4s" %(z+1)
           if supercell in line:
              # Go through all lines for a supercell
              for n in range(ions):
                 for i in range(3): 

                    # CASTEP data is contained into 6 columns with ceil(atoms/2) rows.
                    # line_data stores this data into a single row.
                    line_data = []
                    for f in range(l_castep):                                             
                        new_line = index + (n*3)*l_castep + i*l_castep + f + 1                 
                        line_split = castep_data[new_line].split()
                        #print 'line', new_line , line_split

                        if len(line_split)==8:
                           line_data+=(float(line_split[j]) for j in xrange(2, 8))
                        elif len(line_split)<8 and len(line_split)>3:
                           line_data+=(float(line_split[j]) for j in xrange(0, 6))
                        elif len(line_split)<4:
                           line_data+=(float(line_split[j]) for j in xrange(0, 3))
                     
                    supercell_3nx3n_matrix += [line_data]

                    for m in range(ions):
                        ifc_data[supercell][ions*n+m][i]=[float(supercell_3nx3n_matrix[3*n+i][j]) for j in xrange(m*3, 3*m+3)]

                  
    #print 'cell 1 atom 1',ifc_data['Supercell   1'][1] 



    """
    All of the necessary information is loaded into the dictionary ifc_data. 
    It would be used from now on to write out all permutations.
    """    
    
# Required format:
# line 1: nx*xy*nz*ions
# line 2: atom1 atom2
# lines 3,4,5: IFC matrix for atom 1 and 2
#
# Permutaions use the following pattern: Pick atom 1 and write out its interations with all other atoms
# For example, in a 5x5x5 supercell with 2 atoms, say Na and Cl, in the unit cell, the interation between atoms 1 and 2
# is the interaction between Na in cell (1,1,1) and Na in (2,1,1). 1-125 corresponds to Na(1,1,1) and Na(5,5,5).
# Therefore the interaction between Na and Cl in (1,1,1) is given by the pair 1-126.   
# Supercell is constructured by increasing x first, then y and then z.
      
    permutations = np.empty((nx,ny,nz))
    #print permutations[1,1,1] 
    for m in range(len(ifc_data['Supercell   1'])):
      a = np.append(permutations[0][0][0], [ifc_data['Supercell   1'][m]])
      #permutations[0,0,0]=[ifc_data['Supercell   1'][m]]
      #print m#, ifc_data['Supercell   1'][m]
    #print permutations, 'ab', a    
    


    j=0
    j2=0 
    indeces_array = []
    for every_atom in xrange(0,ions*nx*ny*nz):
        print "ATOM %s" %(every_atom+1)
        j+=1
        j2+=1

        if j2 > nx*ny*nz:
          j2=1
        
        lol, wow = divmod(every_atom, nx)
        lol1, wow1 = divmod(every_atom, nx*ny)
    	for every_ion in xrange(0,ions):
	
        	for z in xrange(0,nz):
 			
                        l2 = z + lol1 + every_ion*nz
                        oz, uz = divmod(l2, nz)
			for y in xrange(0,ny):	
				
                                l1 = y + lol + z*ny + every_ion*nz*ny
                                oy,uy = divmod(l1, ny)
                                for x in xrange(0,nx):
				     
                                     l = (x+j2) + y*nx + z*ny*nx + every_ion*nz*ny*nx
                                     ox,ux = divmod(l-1,nx)
                                                                          
                                     #print j,  (ux+1+nx*uy+uz*ny*nx + every_ion*nz*ny*nx)
                                     indeces_array.append([int(j), int((ux+1+nx*uy+uz*ny*nx + every_ion*nz*ny*nx))])
                                 
                                   
 
    print len(ifc_data) # < gives the number of supercells
    print indeces_array

    
    keylist = ifc_data.keys()
    keylist.sort()
    

    ifc_for_sheng=dict()
    i=0
    # go through all columns
    for x in xrange(0,ions):
        
        # do the process for all copies of atom (x) 
        # [just keep the same column and print its data again]
        for lp in xrange(0,(nx*ny*nz)):
             
	     print 'ATOM %s' % (ions*x + lp + 1) # debug
	     
             # extract data from a column (given by x) for the interaction between
             # atom N with a given atom (defined by x above)
	     for y in xrange(0,ions):
                   # go through all supercells
   		   for key in keylist: 
    			    i+=1       
     			    #print (x*nx+1),i
                            print indeces_array[i-1]
                            # here [1] and [0] define the direction of the force
                            ion_couple = '%4s'%indeces_array[i-1][1] + '%4s' %indeces_array[i-1][0]
                            ifc_for_sheng[ion_couple]= [ifc_data[key][ions*y+x]]
                            
                            print ion_couple
                            
                           #print 'i:', i, 'y:' ,y, 'lp:', lp, 'x:', x
        		    print ifc_data[key][ions*y+x]  
    


    f = open('IFC2nd', 'w')
    f_struct = '%4s' %str(nx*ny*nz*ions) + '\n'
    for key1 in sorted(ifc_for_sheng):
       f_struct += str(key1) + '\n'
       for m in range(3):
         for l in range (3):
         #number2 = '%15s'%ifc_for_sheng[key1][0][m][0]
         #print number2
           f_struct += '%22s' %str("{0:.15f}".format(ifc_for_sheng[key1][0][m][l])) 
         f_struct += '\n'
      # print key1
      # print ifc_for_sheng[key1]
    #print ifc_for_sheng['4   4']
    f.write(f_struct)
    f.close()                         

if __name__ == "__main__":
    import sys
    sys.exit(main())


