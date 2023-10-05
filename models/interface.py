
from . import carbonlab
from ..aux   import convert_xyz2abc, check_file ### for model, 1elec
from ..atoms   import *
from ..aux.env import *
import sys, os, importlib 

def interface():
    return 0

def twoInterfaces(dcenter, lphase, rphase, junc_dist):
    '''
    Make cell with two interfaces
        Receive center phase with 'name', dcenter['name']
        Receive size of center with key 'size', dcenter['center']
        combine lphase + modelling_center + rphase

    dcenter: dictionary for scattering region
        name: key words for scattering region
        size: the size of   scattering region
    '''
    cwd = os.getcwd()
         
    if dcenter['name'] == 'grp':
        center = carbonlab.grp_rect(2,dcenter['size'])
        center.select_all()
        center.rotate(90, axis_dir=(1,0,0), with_cell=True)
        center.translate(-1.42, 0, 0)
        center.wrap_positions()
        center.select_z(-0.1, 1)
        center.delete()
        ### 2.2 Combine left and right parts of electrode
        center_cell = center.get_cell()
        elecL_cell = lphase.get_cell()
        elecL_cell_abc = convert_xyz2abc(elecL_cell[0], elecL_cell[1], elecL_cell[2])
        center_cell_abc = convert_xyz2abc(center_cell[0], center_cell[1], center_cell[2])
        ratio = elecL_cell_abc[0] / center_cell_abc[0]
        center2 = center.adjust_cell_size(ratio, 7)

        center2.select_all()
        center2.translate(0,5,junc_dist+lphase.get_zmax()-center2.get_zmin())
        new = lphase + center2
        rphase.select_all()
        rphase.translate(0,0,junc_dist+new.get_zmax()-rphase.get_zmin())
        new = new + rphase
        elecL_cell[2][2] = new.get_zmax() - new.get_zmin()
        new.set_cell(elecL_cell)
        new.set_vacuum(2.306)

        new_cell = new.get_cell()
        new_zavg = (new.get_zmin()+new.get_zmax())/2
        cell_avg = convert_xyz2abc(new_cell[0], new_cell[1], new_cell[2])[2]/2
        new.select_all()
        new.translate(0,0,cell_avg - new_zavg)
        new.select_all()
        new.sort(option = 'z'); new.set_serials(1)
        
    ### Further coding for Other lab
    
    return new
