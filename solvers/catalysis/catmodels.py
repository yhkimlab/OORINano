
import os, sys, glob, math
from ...atoms import *

'''
only class Catmodels
'''


class Catmodels:
    '''
    atoms   AtomsSystems instance
    '''

    def __init__(self, atoms):
        self.atoms = atoms

    def get_zmax_index(self):
        atoms = self.atoms
        z_axis = []
        for atom in atoms:
            x, y, z = atom.get_position()
            z_axis.append(z)
        maxindex = z_axis.index(max(z_axis))
        return maxindex

    def HER_transition_gen(self, pivot=None, zdist=1.5):
        """
        specify hydrogen atomic position on the catalyst

        Parameters
        ----------
        pivot = list of position or int(atom index)
        zdist   = 1.5
            assign the position [A] of adsorbed hydrogen [x, y, z]
            Ex. pivot = [1.2, 1.35, 23.75]
            if pivot is None, the atomic position is automatically assigned
            along to highest z-axis position (ontop site) with specific distance (dist=1.5)
        """
        atoms = self.atoms
        vHdist     = Vector(  0.,      0.,     zdist)
        if pivot is None: # do not obtain pivot after geometry is disturbed
           top_index = self.get_zmax_index()
           pivot_position = atoms[top_index][0], atoms[top_index][1], atoms[top_index][2]
        else:
           if type(pivot) is int:
                pivot_position = atoms[pivot-1].get_position()
           else:
                pivot_position = (pivot[0], pivot[1], pivot[2])

        vpivot = Vector(pivot_position)
        atomsH = atoms + Atom('H', vpivot + vHdist)
        return atomsH

    def four_electron_transition_gen(self, pivot=None, zdist=1.5, mode='ORR'):
        """
        specify initial atomic position on the catalyst

        Parameters
        ----------
        pivot = int (pivot id), list
            same as HER_transition_gen
            assign the position [A] of adsorbed hydrogen [x, y, z]
            Ex. pivot = [1.2, 1.35, 23.75]
            if pivot is None, the atomic position is automatically assigned
            along to highest z-axis position (ontop site) with specific distance (dist=1.5)
        dist   = 1.5
        vdist   vector distance
        mode = 'ORR' or 'OER'
        
        return 4 atoms images
        """
        atoms = self.atoms
        vO1dist     = Vector(  0.,      0.,     zdist)
        vO2dist     = Vector(-1.0,      0.4,    zdist + 0.6 )
        vO2Hdist    = Vector(-0.403,    1.054,  zdist + 0.733)
        vHdist      = Vector(-1.350,    1.196,  zdist + 0.480) 
        vH2dist     = Vector(   0.,     0.,     zdist+0.971)   
        
        if pivot is None: # find center
            top_index = self.get_zmax_index()
            pivot_position = atoms[top_index][0], atoms[top_index][1], atoms[top_index][2]
            
        else:
            if type(pivot) is int:
                pivot_position = atoms[pivot-1].get_position()
                #print(f"pivot_position: {pivot_position}")
            else:
                pivot_position = (pivot[0], pivot[1], pivot[2])
        vpivot = Vector(pivot_position)
        
        atomsO2  = atoms \
                   + Atom('O', vpivot+vO1dist) \
                   + Atom('O', vpivot+vO2dist) 

        atomsOOH = atoms \
                   + Atom('O', vpivot+vO1dist) \
                   + Atom('O', vpivot+vO2Hdist) \
                   + Atom('H', vpivot+vHdist)
        atomsO   = atoms \
                   + Atom('O', vpivot+vO1dist) 

        atomsOH  = atoms \
                   + Atom('O', vpivot+vO1dist) \
                   + Atom('H', vpivot+vH2dist)

        if mode == 'ORR':
            return atomsO2, atomsOOH, atomsO, atomsOH 
        
        elif mode == 'OER':
            return atomsOH, atomsO, atomsOOH
        
    



