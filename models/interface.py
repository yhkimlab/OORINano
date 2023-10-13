
from . import carbonlab
from ..ncio   import convert_xyz2abc ### for model, 1elec
from ..atoms   import *
from ..aux.env import *
import os 

def interface():
    return 0

def twoInterfaces(dict_middle, lphase, rphase, junc_dist):
    '''
    Make cell with two interfaces
        Receive middle phase with 'name', dict_middle['name']
        Receive size of middle with key 'size', dict_middle['middle']
        combine lphase + modelling_center + rphase

    dict_middle: dictionary for scattering region
        name: key words for scattering region
        size: the size of   scattering region
    '''
    cwd = os.getcwd()
         
    if dict_middle['name'] == 'grp':
        middle = carbonlab.grp_rect(2,dict_middle['size'])
        middle.select_all()
        middle.rotate(90, axis_dir=(1,0,0), with_cell=True)
        middle.translate(-1.42, 0, 0)
        middle.wrap_positions()
        middle.select_z(-0.1, 1)
        middle.delete()
        
        middle_cell = middle.get_cell()
        elecL_cell = lphase.get_cell()
        elecL_cell_abc = convert_xyz2abc(elecL_cell[0], elecL_cell[1], elecL_cell[2])
        middle_cell_abc = convert_xyz2abc(middle_cell[0], middle_cell[1], middle_cell[2])
        ratio = elecL_cell_abc[0] / middle_cell_abc[0]
        middle2 = middle.adjust_cell_size(ratio, 7)

        middle2.select_all()
        middle2.translate(0,5,junc_dist+lphase.get_zmax()-middle2.get_zmin())
        junc_model = lphase + middle2
        rphase.select_all()
        rphase.translate(0,0,junc_dist+junc_model.get_zmax()-rphase.get_zmin())
        junc_model = junc_model + rphase
        elecL_cell[2][2] = junc_model.get_zmax() - junc_model.get_zmin()
        junc_model.set_cell(elecL_cell)
        junc_model.set_vacuum(2.306)

        junc_cell = junc_model.get_cell()
        junc_zavg = (junc_model.get_zmin()+junc_model.get_zmax())/2
        cell_avg = convert_xyz2abc(junc_cell[0], junc_cell[1], junc_cell[2])[2]/2
        junc_model.select_all()
        junc_model.translate(0,0,cell_avg - junc_zavg)
        junc_model.select_all()
        junc_model.sort(option = 'z'); junc_model.set_serials(1)
        
    ### Further coding for Other lab
    
    return junc_model
