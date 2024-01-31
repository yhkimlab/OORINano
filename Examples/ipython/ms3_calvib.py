# coding: utf-8
from oorinano.calculator.vasp import readAtomicStructure as read_geo
from oorinano.calculator.vasp import Vasp
atoms = read_geo("CONTCAR_test_catH")
calc = Vasp(atoms)
zpe, ts = calc.get_vibration_energy(output_name='OUTCAR_test_catH_vib')
zpe
ts
calc.get_vibration_spectrum(output_name='OUTCAR_test_catH_vib')
