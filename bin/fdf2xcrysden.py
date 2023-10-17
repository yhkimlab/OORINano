#!/usr/bin/env python 
from oorinano import *
from oorinano.simulators import siesta as s2
from oorinano.aux import vis

def fdf2xcrysden(fdf_name):
    at = s2.readAtomicStructure(fdf_name)
    vis.show_xcrysden(at)

if __name__ == '__main__':

    import sys
    
    usage = """
    Usage :
      xyz2xcrysden.py <xyz file>
    """ 

    if len(sys.argv) != 2:
        print(__doc__); print(usage)
        sys.exit()
    else:
        fdf_name = sys.argv[1]

    fdf2xcrysden(fdf_name)
