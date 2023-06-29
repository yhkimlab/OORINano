import os, sys
import copy
import numpy as np
import nanocore as nc
from nanocore import siesta as s2
from time import sleep

from nanocore.atoms import AtomsSystem, Vector

option = 'A'  # ABABA
delta = [2.88, 3.269] # interlayer, G-TMDC
gate = False

def parse_option():
    temp = ' '.join(option)
    res = temp.split()
    return res

def change_model(model:AtomsSystem, opt):
    new = []
    if opt == 'A':
        symbol = list(set(model.get_symbols()))
        assert len(symbol) == 2
        for atom in model:
            if atom.get_symbol() == symbol[0]:
                x1 = atom.get_position()
            elif atom.get_symbol() == symbol[1]:
                x2 = atom.get_position()
        delta = x1 - x2
        delta.x3 = 0
        new = model
        new.select_all()
        new.translate(delta.x1, delta.x2, delta.x3)
    elif opt == 'B':
        symbol = list(set(model.get_symbols()))
        assert len(symbol) == 2
        for atom in model:
            if atom.get_symbol() == symbol[0]:
                x1 = atom.get_position()
            elif atom.get_symbol() == symbol[1]:
                x2 = atom.get_position()
        for atom in model:
            if atom.get_symbol() == symbol[0]:
                vector = atom.get_position()
                vector.x1 = x2.x1
                vector.x2 = x2.x2
                atom.set_position(vector)
                new.append(atom)
            elif atom.get_symbol() == symbol[1]:
                vector = atom.get_position()
                vector.x1 = x1.x1
                vector.x2 = x1.x2
                atom.set_position(vector)
                new.append(atom)
        new = AtomsSystem(new)
        new.set_cell(model.get_cell())
    return new
        

if __name__ =='__main__':
    model = s2.read_fdf("MoS2.fdf")
    elec = s2.read_fdf("GRP.fdf")
    opt = parse_option()
    
    m1, m2 = elec.get_cell_match(model, 'x', 'x', 5)
    modelB = change_model(model, 'B')
    elec2 = elec * (m1, m1, 1)
    modelA = model * (m2, m2, 1)
    modelB = modelB * (m2, m2, 1)

    ec = elec2.get_cell()
    mc = modelA.get_cell()
    elec_cell = nc.io.convert_xyz2abc(ec[0], ec[1], ec[2])
    model_cell = nc.io.convert_xyz2abc(mc[0], mc[1], mc[2])
    ratio = elec_cell[0] / model_cell[0]
    modelA2 = modelA.adjust_cell_size(ratio, 4)
    modelB2 = modelB.adjust_cell_size(ratio, 4)

    new = copy.deepcopy(elec2)
    for i, op in enumerate(opt):
        if op == 'A':
            model2 = modelA2
        elif op == 'B':
            model2 = modelB2
        model_zmin = model2.get_zmin()
        new_zmax = new.get_zmax()
        model2.select_all()
        if i == 0:
            model2.translate(0,0,delta[1]+new_zmax-model_zmin)
            new = new + model2
        else:
            model2.translate(0,0,delta[0]+new_zmax-model_zmin)
            new = new + model2

    new_zmax = new.get_zmax()
    elec2_zmin = elec2.get_zmin()
    elec2.select_all()
    elec2.translate(0,0,delta[1]+new_zmax-elec2_zmin)
    new = new + elec2

    if gate == True:
        Au = nc.atoms.Atom('Au', Vector(0,0,new.get_zmax() + 11.53))
        Au_slab = AtomsSystem([Au], elec.get_cell())
        Au_slab = Au_slab * (m1, m1, 1)
        new = new + Au_slab

    new.set_cell(modelA2.get_cell())
    new_zavg = (new.get_zmin()+new.get_zmax())/2
    cell_avg = model_cell[2]/2
    new.select_all()
    new.translate(0,0,cell_avg - new_zavg)

    sim = s2.Siesta(new)
    sim.write_struct()

