import nanocore as nc
from nanocore.models import carbonlab
from nanocore.simulator.siesta import Siesta
from nanocore.io import convert_xyz2abc, write_struct, read_struct

elecL = read_struct("initial/elecL.fdf")
elecR = read_struct("initial/elecR.fdf")

max_unit = 12
min_unit = 3
junction_distance = 1.96 # given parameter
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
    new.set_vacuum(20.0)

    new_cell = new.get_cell()
    new_zavg = (new.get_zmin()+new.get_zmax())/2
    cell_avg = convert_xyz2abc(new_cell[0], new_cell[1], new_cell[2])[2]/2
    new.select_all()
    new.translate(0,0,cell_avg - new_zavg)
    new.select_all()
    new.sort(option = 'z'); new.set_serials(1)
    write_struct(new, fname = f"cnt_{i}.fdf")
    # write_poscar(new, f"cnt_{i}.poscar")
    
