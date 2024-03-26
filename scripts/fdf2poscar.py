#!/usr/bin/env python
from oorinano import *
from oorinano.calculator.siesta import readAtomicStructure  as read_geo
from oorinano.calculator.vasp   import writeAtomicStructure as write_geo
import sys

fname = sys.argv[1]
product = fname.split('.')[0]
atom = read_geo(fname)

write_geo(atom, f'{product}.poscar')


