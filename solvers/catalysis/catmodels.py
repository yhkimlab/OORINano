
import os, sys, glob, math
from ...atoms import *

'''
only class Catmodeling
'''


class Catmodeling:

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

    def HER_transition_gen(self, active=None, dist=1.5):
        """
        specify hydrogen atomic position on the catalyst

        Parameters
        ----------
        active = list
        dist   = 1.5
            assign the position [A] of adsorbed hydrogen [x, y, z]
            Ex. active = [1.2, 1.35, 23.75]
            if active is None, the atomic position is automatically assigned
            along to highest z-axis position (ontop site) with specific distance (dist=1.5)
        """
        atoms = self.atoms

        if active is None:
           top_index = self.get_zmax_index()
           top_x, top_y, top_z = atoms[top_index][0], atoms[top_index][1], atoms[top_index][2]
           atomsH = atoms + Atom('H', (top_x, top_y, top_z+dist))
        else:
           H_x, H_y, H_z = active[0], active[1], active[2]
           atomsH = atoms + Atom('H', (H_x, H_y, H_z))

        return atomsH

    def four_electron_transition_gen(self, active=None, dist=1.5, mode='ORR'):
        """
        specify initial atomic position on the catalyst

        Parameters
        ----------
        active = list
        dist   = 1.5
            assign the position [A] of adsorbed hydrogen [x, y, z]
            Ex. active = [1.2, 1.35, 23.75]
            if active is None, the atomic position is automatically assigned
            along to highest z-axis position (ontop site) with specific distance (dist=1.5)
        
        mode = 'ORR' or 'OER'
        active = list
            same as HER_transition_gen
        return 4 atoms images
        """
        atoms = self.atoms
        
        if active is None:
           top_index = self.get_zmax_index()
           top_x, top_y, top_z = atoms[top_index][0], atoms[top_index][1], atoms[top_index][2]
           active_position = (top_x, top_y, top_z+dist)
        else:
           active_position = (active[0], active[1], active[2])

        atomsO2  = atoms \
                   + Atom('O', (active_position[0], active_position[1], active_position[2])) \
                   + Atom('O', (active_position[0]-1.000, active_position[1]+0.400, active_position[2]+0.600)) 

        atomsOOH = atoms \
                   + Atom('O', (active_position[0], active_position[1], active_position[2])) \
                   + Atom('O', (active_position[0]-0.403, active_position[1]+1.054, active_position[2]+0.733)) \
                   + Atom('H', (active_position[0]-1.350, active_position[1]+1.196, active_position[2]+0.480))
        atomsO   = atoms \
                   + Atom('O', (active_position[0], active_position[1], active_position[2])) 

        atomsOH  = atoms \
                   + Atom('O', (active_position[0], active_position[1], active_position[2])) \
                   + Atom('H', (active_position[0], active_position[1], active_position[2]+0.971))

        if mode == 'ORR':
            return atomsO2, atomsOOH, atomsO, atomsOH 
        
        elif mode == 'OER':
            return atomsOH, atomsO, atomsOOH
        
    



