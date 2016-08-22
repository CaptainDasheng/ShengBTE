#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage ./castep2shengbte.py <seedname>
# where seedname is the name of the .castep file.
#
#
# v0.2 11 May 16 - now loads castep IFC data into a dictionary
# for each supercell. Each supercell contains ions^2 IFC matrices.
#   
# 
# v0.3 08 Aug 16 - a working prototype, which generates the 2nd order IFC file
# TODO: Now it is time to add lines which create the control file
#
#
# v1.0 22 Aug 16 - the program creates the 2nd order IFC and CONTROL files.
# CONTROL file is missing epsilon_0 - the dielectric tensor information.


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
    print '||                     version 1.0                     ||'
    print '||                     22 Aug 2016                     ||'
    print '||-----------------------------------------------------||'
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
       # Avoid ugly errors
       print '||        Usage: castep2shengbte.py <seedname>         ||'
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
            print '||        Usage: castep2shengbte.py <seedname>         ||'
            print '========================================================='
            sys.exit()


    """
    Return all the relevant information contained in a .castep file.
    """
    castep_file = prefix + '.castep'
    castep_file = open(castep_file, 'r')
    castep_data = castep_file.readlines()
    castep_file.close()
   
    for line in castep_data:
      if 'Total number of ions in cell' in line:
         ions = int(line.split()[7])
      elif 'Total number of species in cell' in line:
         species = int(line.split()[7])
      elif 'MP grid size for SCF calculation is' in line:  
         nx = int(line.split()[7])
         ny = int(line.split()[8])
         nz = int(line.split()[9])
         supercell_size = nx*ny*nz

    # Get unit cell from <seedname>.castep.
    lat_vectors = []
    for index, line in enumerate(castep_data):
         # Crystal lattice in (A) units. <seedname>.bands file also contains
         # this information. However, it uses Bohr units and this creates 
         # some problems when symmetry operations are generated with spglib. 
         if 'Real Lattice(A)' in line:  
                           start = index + 1
                           for j in range (3):
                             lat_vectors.append(
                                              [float(castep_data[start].split()[0]), 
                                               float(castep_data[start].split()[1]),
                                               float(castep_data[start].split()[2])])
                             start += 1  
                           break # avoid double counting    

    # Species names
    species_names = []  
    for index, line in enumerate(castep_data): 
	 if ' Mass of species in AMU' in line:
			   start1 = index + 1
			   for j1 in range (int(species)):
                               species_names.append(str(castep_data[start1].split()[0]))
                               start1 += 1
                           break # avoid double counting

    # Species coordinates
    positions = []
    for index, line in enumerate(castep_data):
         if 'Cell Contents' in line:
            for i in range(0, ions):
                positions.append([str(castep_data[index+10+i].split()[1]),
				  float(castep_data[index+10+i].split()[3]), 
                                  float(castep_data[index+10+i].split()[4]), 
                                  float(castep_data[index+10+i].split()[5])])
            break # avoid double counting   

    # Create a list for types which is used in the CONTROL file
    types = []
    for number_ions in range(len(positions)):
		for position, item in enumerate(species_names):
                  if positions[number_ions][0] in item:
                    types.append(position+1)

    # Born Effective Charges
    born_charges = []
    for index, line in enumerate(castep_data):
         if 'Born Effective Charges' in line:
            for i in range(0, ions):
              for j in range(3): 
                if len(castep_data[index+2+i*3+j].split()) >= 4:
                    born_charges.append([float(castep_data[index+2+i*3+j].split()[2]), 
                                  float(castep_data[index+2+i*3+j].split()[3]), 
                                  float(castep_data[index+2+i*3+j].split()[4])])
                   
                elif len(castep_data[index+2+i*3+j].split())<=3:
                  born_charges.append([float(castep_data[index+2+i*3+j].split()[0]), 
                                  float(castep_data[index+2+i*3+j].split()[1]), 
                                  float(castep_data[index+2+i*3+j].split()[2])])
                   
              break # avoid double counting   
                    
    ifc_data=dict()   
    for i in range(supercell_size):
        supercell = "Supercell" + "%4s" %(i+1)
        ifc_data[supercell]=[]         
        for n in range(ions*ions):
           ifc_data[supercell].append(np.empty((3,3)))
        
    # Define a number which gives the number of lines we have per atom per space index in the force constants matrix
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


    """
    All of the necessary information is loaded into the dictionary ifc_data. 
    It would be used from now on to write out all permutations.
    """    
    
# Required format in the FORCE_CONSTANTS_2ND file:
# line 1: nx*xy*nz*ions
# line 2: atom1 atom2
# lines 3,4,5: IFC matrix for atom 1 and 2
#
# Permutaions use the following pattern: Pick atom 1 and write out its interations with all other atoms
# For example, in a 5x5x5 supercell with 2 atoms, e.g. NaCl in the unit cell, the interaction between atoms 1 and 2
# is the interaction between Na in cell (1,1,1) and Na in (2,1,1). 1-125 corresponds to Na(1,1,1) - Na(5,5,5), respectively.
# Therefore the interaction between Na and Cl in (1,1,1) is given by the pair 1-126.   
# The supercell is constructured by increasing x first, then y and then z.

    ion_reset=0 
    indeces_array = []
    for every_atom in xrange(0,ions*nx*ny*nz):
        # pick 1 ion, go through all cells and then
        # reset the number to 1 for the next ion
        ion_reset+=1
        if ion_reset > nx*ny*nz:
          ion_reset=1
        
        divnx, dummy1 = divmod(every_atom, nx)
        divnxny, dummy2 = divmod(every_atom, nx*ny)

    	for every_ion in xrange(0,ions):	
        	for z in xrange(0,nz):
 			counter_z = z + divnxny + every_ion*nz
                        dummy3, modnz = divmod(counter_z, nz)
			
			for y in xrange(0,ny):					
                                counter_y = y + divnx + z*ny + every_ion*nz*ny
                                dummy4, modny = divmod(counter_y, ny)
                                
				for x in xrange(0,nx):				     
                                     counter_x = (x+ion_reset) + y*nx + z*ny*nx + every_ion*nz*ny*nx
                                     dummy5, modnx = divmod(counter_x-1,nx)                                     
                                     indeces_array.append([int(every_atom+1), int((modnx+1+nx*modny+modnz*ny*nx + every_ion*nz*ny*nx))])
                                 
                                    
    keylist = ifc_data.keys()
    keylist.sort()
    ifc_for_sheng=dict()
    i=0
    # go through all columns
    for x in xrange(0,ions):
        
        # do the process for all copies of atom (x) 
        # [just keep the same column and print its data again]
        for lp in xrange(0,(nx*ny*nz)):
             
             # extract data from a column (given by x) for the interaction between
             # atom N with a given atom (defined by x above)
	     for y in xrange(0,ions):
                   # go through all supercells
   		   for key in keylist: 
    			    i+=1       
     			    # here [i-1][0] and [i-1][1] define the direction of the forces
                            ion_couple = '%4s' %indeces_array[i-1][0] + '%4s' %indeces_array[i-1][1]
                            ifc_for_sheng[ion_couple]= [ifc_data[key][ions*y+x]]
                            
    # Create the 2nd order file
    f = open('FORCE_CONSTANTS_2ND', 'w') # use test name 'IFC2nd' for now
    f_struct = '%4s' %str(nx*ny*nz*ions) + '\n'
    for key1 in sorted(ifc_for_sheng):
       f_struct += str(key1) + '\n'
       for m in range(3):
         for l in range (3):
           f_struct += '%22s' %str("{0:.15f}".format(ifc_for_sheng[key1][0][m][l])) 
         f_struct += '\n'
    f.write(f_struct)
    f.close()   

#################################################################################################################################################

    # Create the control file   
    f = open('CONTROL', 'w') # use test name 'CONTROL-CASTEP' for now

    f_control  = '&allocations' +'\n'
    f_control += '\tnelements=' + '%s,' %str(species) + '\n'  
    f_control += '\tnatoms=' + '%s,' %str(ions) + '\n'
    f_control += '\tngrid(:)=' + '%s ' %str(nx) + '%s ' %str(ny) +'%s' %str(nz) +'\n'
    # f_control += '\tnorientations=' + '%s' %str() + '\n' # (integer, default=0): number of orientations along which to study nanowires (not needed in CASTEP for now)
    f_control += '&end' +'\n'

    f_control += '&crystal' +'\n'
    # f_control += '\tlfactor=' + '%s,' %str(species) + '\n' # not needed in CASTEP, it is possible to use latt vecs straight away
    f_control += '\tlattvec(:,1)=' + '%s ' %str(lat_vectors[0][0]/10) + '%s ' %str(lat_vectors[0][1]/10)+ '%s,' %str(lat_vectors[0][2]/10) + '\n' 
    f_control += '\tlattvec(:,2)=' + '%s ' %str(lat_vectors[1][0]/10) + '%s ' %str(lat_vectors[1][1]/10)+ '%s,' %str(lat_vectors[1][2]/10) + '\n' 
    f_control += '\tlattvec(:,3)=' + '%s ' %str(lat_vectors[2][0]/10) + '%s ' %str(lat_vectors[2][1]/10)+ '%s,' %str(lat_vectors[2][2]/10) + '\n' 
    f_control += '\telements='
    for el in range(len(species_names)):
         f_control +='"%s" ' %str(species_names[el])
    f_control += '\n'
    f_control += '\ttypes='
    for el_num in range(len(types)):
		f_control += ' %s' %str(types[el_num])
    f_control += ',\n'  
    
    for atoms in range(len(positions)):
       		f_control += '\tpositions(:,%s)=' %str(atoms+1) 
                for i in range(1,4):
                   f_control += ' %s' %str("{0:.8f}".format(positions[atoms][i]))
                f_control += ',\n'  

    for atoms in range(len(positions)):
          for i in range(1,4):
       		f_control += '\tborn(:,%s,' %str(i) + '%s)=' %str(atoms+1)
                for j in range(3):
                   f_control += ' %s' %str("{0:.5f}".format(born_charges[atoms*3 + (i-1)][j]))
                f_control += ',\n'  

    f_control += '\tscell(:)=' + '%s ' %str(nx) + '%s ' %str(ny) +'%s' %str(nz) +'\n'
    f_control += '&end' +'\n'


    f_control += '&parameters' +'\n'
    f_control += '\tT=300.'+'\n'
    f_control += '\tscalebroad=1.0'+'\n'
    f_control += '&end' +'\n'

    f_control += '&flags' +'\n'
    f_control += '\tnonanalytic=.TRUE.'+'\n'

    f_control += '&end' +'\n'

    f.write(f_control)
    f.close()                      
####################################################################################################################################################


if __name__ == "__main__":
    import sys
    sys.exit(main())


