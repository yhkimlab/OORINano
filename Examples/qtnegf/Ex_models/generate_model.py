from oorinano.models import carbonlab
from oorinano.calculator.siesta import Siesta, readAtomicStructure, writeAtomicStructure
from oorinano.calculator.siesta import default_location
from oorinano.utils.auxil import convert_xyz2abc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-M", "--max", dest="max_unit", action="store", type=int, default=12)
parser.add_argument("-m", "--min", dest="min_unit", action="store", type=int, default=3)
parser.add_argument("-d", "--dist", dest="dist", action="store", type=float, default=1.96)
parser.add_argument("-l", "--left", dest="left", action="store")
parser.add_argument("-r", "--right", dest="right", action="store" )
args = parser.parse_args()

if args.left:
    path_elecL = args.left
    path_elecR = args.right
else:    
    path_elecL = default_location + '/model/channel/elecL.fdf'
    path_elecR = default_location + '/model/channel/elecR.fdf'

elecL = readAtomicStructure(path_elecL)
elecR = readAtomicStructure(path_elecR)

max_unit = args.max_unit
min_unit = args.min_unit
junction_distance = args.dist # given parameter
for i in range(min_unit, max_unit, 1):
    cnt = carbonlab.grp_rect(2,i)
    cnt.select_all()
    cnt.rotate(90, axis_dir=(1,0,0), with_cell=True)
    cnt.translate(-1.42, 0, 0)
    cnt.wrap_positions()
    cnt.select_z(-0.1, 1)
    cnt.delete()

    cnt_cell = cnt.get_cell()
    elecL_cell = elecL.get_cell()
    elecL_cell_abc = convert_xyz2abc(elecL_cell[0], elecL_cell[1], elecL_cell[2])
    cnt_cell_abc = convert_xyz2abc(cnt_cell[0], cnt_cell[1], cnt_cell[2])
    ratio = elecL_cell_abc[0] / cnt_cell_abc[0]
    cnt2 = cnt.adjust_cell_size(ratio, 7)

    cnt2.select_all()
    cnt2.translate(0,5,junction_distance+elecL.get_zmax()-cnt2.get_zmin())
    new = elecL + cnt2
    elecR.select_all()
    elecR.translate(0,0,junction_distance+new.get_zmax()-elecR.get_zmin())
    new = new + elecR
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
    writeAtomicStructure(new, fname = f"cnt_{i}.fdf")
    # write_poscar(new, f"cnt_{i}.poscar")
    
