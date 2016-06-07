#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  thirdorder, help compute anharmonic IFCs from minimal sets of displacements
#  Copyright (C) 2012-2014 Wu Li <wu.li.phys2011@gmail.com>
#  Copyright (C) 2012-2014 Jesús Carrete Montaña <jcarrete@gmail.com>
#  Copyright (C) 2012-2014 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
#  Copyright (C) 2014      Antti J. Karttunen <antti.j.karttunen@iki.fi>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.




#  Using thirdorder_vasp.py as a basis to create thirdorder_castep.py
#  1) start by converting [read_POSCAR] to read_CASTEP
#  2) make the [write_POSCAR] output to a suitable castep form
#  3) Find out what needs to be converted for the reap command 

import os.path
import glob
try:
    from lxml import etree as ElementTree
    xmllib="lxml.etree"
except ImportError:
    try:
        import xml.etree.cElementTree as ElementTree
        xmllib="cElementTree"
    except ImportError:
        import xml.etree.ElementTree as ElementTree
        xmllib="ElementTree"
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO
try:
    import hashlib
    hashes=True
except ImportError:
    hashes=False

import sys
sys.path.insert(0, '../')

import thirdorder_core
from thirdorder_common import *

# Convert this to castep. Will use a cell file for the moment

########################################################################
# CASTEP

def read_POSCAR(directory):
    """
    Return all the relevant information contained in a .cell file.
    """
    with dir_context(directory):
        nruter=dict()
        nruter["lattvec"]=np.empty((3,3))
        nruter["elements"]=[]
        
        f=open("castep.cell","r")
        castep_cell = f.readlines()
        

        atoms_list = []
        for index, line in enumerate(castep_cell):
           if '%BLOCK lattice_cart' in line:              
              for i in xrange(3):
                  new_line = index + i + 1
                  nruter["lattvec"][:,i]=[float(j) for j in castep_cell[new_line].split()]
           #elif '%BLOCK lattice_abc' in line: # add this part later
           elif '%BLOCK positions_frac' in line:
                index_start = index
                print 'start block', index_start
           elif '%ENDBLOCK positions_frac' in line:
                index_end = index
                print 'end index', index_end
        for i in range(index_start+1, index_end):
               atoms_list.append(castep_cell[i].split())
        atoms_list = filter(None,atoms_list)
        atoms_list.sort(key=lambda tup: tup[0])  # for now
        print ' TEST' ,atoms_list
        natoms1 = len(atoms_list)
        nruter["positions"]=np.empty((3,natoms1))
        for i in range(natoms1):
             
            nruter["positions"][:,i]=[float(atoms_list[i][j]) for j in xrange(1,4)]
            nruter["elements"].append(str(atoms_list[i][0]))
        print 'nruter["elements"]', nruter["elements"]
        nruter["elements"].sort(key=lambda tup: tup[0])
        test = nruter["elements"]
        nruter["elements"]=list(set(nruter["elements"]))
        nruter["elements"].sort()
        print 'nruter["elements"]', nruter["elements"]
        

     
        nruter["numbers"]= np.array([int(test.count(nruter["elements"][i])) for i in range(len(nruter["elements"]))],dtype=np.intc)  
        nruter["types"]=[]
        print 'nruter["numbers"]', nruter["numbers"]
        for i in xrange(len(nruter["numbers"])):
             nruter["types"]+=[i]*nruter["numbers"][i]
       
      #  print 'test2',  nruter1["elements"], nruter1["numbers"], test.count(nruter1["elements"][1])    
        
    return nruter

def write_POSCAR(poscar,filename):
    """
    Write the contents of poscar to filename.
    """
    global hashes
    f=StringIO.StringIO()
    #f.write("1.0\n")
    f.write("%BLOCK lattice_cart\n")
    for i in xrange(3):
        f.write("{0[0]:>20.15f} {0[1]:>20.15f} {0[2]:>20.15f}\n".format(
            (poscar["lattvec"][:,i]).tolist()))
    f.write("%ENDBLOCK lattice_cart\n")
    f.write("\n")
    f.write("%BLOCK positions_frac\n")
    
    k = 0
    for i in xrange(len(poscar["numbers"])):
        
        for j in xrange(poscar["numbers"][i]):
           l = k + j
           f.write("{0}".format("".join(poscar["elements"][i]))) 

           f.write("{0[0]:>20.15f} {0[1]:>20.15f} {0[2]:>20.15f}\n".format(
            poscar["positions"][:,l].tolist()))
        k += j + 1
    f.write("%ENDBLOCK positions_frac\n")

    with open(filename,"w") as finalf:
       
        finalf.write(f.getvalue())
    f.close()



def normalize_SPOSCAR(sposcar):
    """
    Rearrange sposcar, as generated by gen_SPOSCAR, so that it is in
    valid VASP order, and return the result.
    """
    nruter=copy.deepcopy(sposcar)
    # Order used internally (from most to least significant):
    # k,j,i,iat For VASP, iat must be the most significant index,
    # i.e., atoms of the same element must go together.
    indices=np.array(xrange(nruter["positions"].shape[1])).reshape(
        (sposcar["nc"],sposcar["nb"],sposcar["na"],-1))
    indices=np.rollaxis(indices,3,0).flatten().tolist()
    nruter["positions"]=nruter["positions"][:,indices]
    nruter["types"].sort()
    return nruter


# CASTEP
def read_forces(filename):
    """
    Read a set of forces on atoms from filename, presumably in
    vasprun.xml format.
    """
   
    f=open(filename,"r")
    castep_forces = f.readlines()
    f.close()    

    nruter = []
    for index, line in enumerate(castep_forces):
         if 'Total number of ions in cell' in line:
            n_atoms = int(line.split()[7])
         if 'Cartesian components (eV/A)' in line:
            starting_line = index + 4
            for i in range(n_atoms):
                
                f = starting_line + i
                nruter.append([float(castep_forces[f].split()[m]) for m in range(3,6)])
    
    nruter=np.array(nruter,dtype=np.double)
    
    return nruter


def build_unpermutation(sposcar):
    """
    Return a list of integers mapping the atoms in the normalized
    version of sposcar to their original indices.
    """
    indices=np.array(xrange(sposcar["positions"].shape[1])).reshape(
        (sposcar["nc"],sposcar["nb"],sposcar["na"],-1))
    indices=np.rollaxis(indices,3,0).flatten()
    
    return indices.argsort().tolist()


if __name__=="__main__":
    if len(sys.argv)!=6 or sys.argv[1] not in ("sow","reap"):
        sys.exit("Usage: {0} sow|reap na nb nc cutoff[nm/-integer]".format(sys.argv[0]))
    action=sys.argv[1]
    na,nb,nc=[int(i) for i in sys.argv[2:5]]
    if min(na,nb,nc)<1:
        sys.exit("Error: na, nb and nc must be positive integers")
    if sys.argv[5][0]=="-":
        try:
            nneigh=-int(sys.argv[5])
        except ValueError:
            sys.exit("Error: invalid cutoff")
        if nneigh==0:
            sys.exit("Error: invalid cutoff")
    else:
        nneigh=None
        try:
            frange=float(sys.argv[5])
        except ValueError:
            sys.exit("Error: invalid cutoff")
        if frange==0.:
            sys.exit("Error: invalid cutoff")
    print "Reading POSCAR"
    poscar=read_POSCAR(".")
    natoms=len(poscar["types"])
    print "Analyzing the symmetries", poscar, natoms
    symops=thirdorder_core.SymmetryOperations(
        poscar["lattvec"],poscar["types"],
        poscar["positions"].T,SYMPREC)
    print "- Symmetry group {0} detected".format(symops.symbol)
    print "- {0} symmetry operations".format(symops.translations.shape[0])
    print "Creating the supercell"
    sposcar=gen_SPOSCAR(poscar,na,nb,nc)
    ntot=natoms*na*nb*nc
    print "Computing all distances in the supercell"
    dmin,nequi,shifts=calc_dists(sposcar)
    if nneigh!=None:
        frange=calc_frange(poscar,sposcar,nneigh,dmin)
        print "- Automatic cutoff: {0} nm".format(frange)
    else:
        print "- User-defined cutoff: {0} nm".format(frange)
    print "Looking for an irreducible set of third-order IFCs"
    wedge=thirdorder_core.Wedge(poscar,sposcar,symops,dmin,
                                nequi,shifts,frange)
    print "- {0} triplet equivalence classes found".format(wedge.nlist)
    list4=wedge.build_list4()
    nirred=len(list4)
    nruns=4*nirred
    print "- {0} DFT runs are needed".format(nruns)
    if action=="sow":
        print sowblock
        print "Writing undisplaced coordinates to 3RD.CASTEP"
        write_POSCAR(normalize_SPOSCAR(sposcar),"3RD.CASTEP.cell")
        width=len(str(4*(len(list4)+1)))
        namepattern="3RD.CASTEP.{{0:0{0}d}}.cell".format(width)
        print "Writing displaced coordinates to 3RD.CASTEP.*"
        for i,e in enumerate(list4):
            for n in xrange(4):
                isign=(-1)**(n//2)
                jsign=-(-1)**(n%2)
                # Start numbering the files at 1 for aesthetic
                # reasons.
                number=nirred*n+i+1
                dsposcar=normalize_SPOSCAR(
                    move_two_atoms(sposcar,
                                   e[1],e[3],isign*H,
                                   e[0],e[2],jsign*H))
                filename=namepattern.format(number)
                write_POSCAR(dsposcar,filename)
    else:
        print reapblock
        print "XML ElementTree implementation: {0}".format(xmllib)
        print "Waiting for a list of vasprun.xml files on stdin"
        filelist=[]
        for l in sys.stdin:
            s=l.strip()
            if len(s)==0:
                continue
            filelist.append(s)
        nfiles=len(filelist)
        print "- {0} filenames read".format(nfiles)
        if nfiles!=nruns:
            print 'sys.exit disabled for a moment'
            # sys.exit("Error: {0} filenames were expected".
           #          format(nruns))
        for i in filelist:
            if not os.path.isfile(i):
                sys.exit("Error: {0} is not a regular file".
                         format(i))
        print "Reading the forces"
        print 'sposcar',sposcar
        p=build_unpermutation(sposcar)
        print 'p',p
        forces=[]
        for i in filelist:
            print read_forces(i)
            forces.append(read_forces(i)[p,:])
            print "- {0} read successfully".format(i)
            res=forces[-1].mean(axis=0)
            print "- \t Average force:"
            print "- \t {0} eV/(A * atom)".format(res)
        print "Computing an irreducible set of anharmonic force constants"
        phipart=np.zeros((3,nirred,ntot))
        for i,e in enumerate(list4):
            for n in xrange(4):
                isign=(-1)**(n//2)
                jsign=-(-1)**(n%2)
                number=nirred*n+i
                phipart[:,i,:]-=isign*jsign*forces[number].T
        phipart/=(400.*H*H)
        print "Reconstructing the full array"
        phifull=thirdorder_core.reconstruct_ifcs(phipart,wedge,list4,poscar,sposcar)
        print "Writing the constants to FORCE_CONSTANTS_3RD"
        write_ifcs(phifull,poscar,sposcar,dmin,nequi,shifts,frange,"FORCE_CONSTANTS_3RD")
    print doneblock
