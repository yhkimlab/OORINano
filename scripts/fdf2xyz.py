#!/usr/bin/env python
from oorinano.calculator.siesta import readAtomicStructure  as read_geo
from oorinano import rw

import sys

fname = sys.argv[1]

atoms = read_geo(fname)
rw.write_xyz('STRUCT.xyz', atoms, comm=None, append=False)

