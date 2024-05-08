#!/usr/bin/env python

from oorinano.calculator.vasp   import readAtomicStructure  as read_geo
from oorinano.calculator.siesta import writeAtomicStructure as write_geo
import sys

fname = sys.argv[1]

atoms = read_geo(fname)

write_geo(atoms)    # fname = 'STRUCT.fdf' in class siesta

