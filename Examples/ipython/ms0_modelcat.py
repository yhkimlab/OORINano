# coding: utf-8
from oorinano import *
at1 = Atom('C', [0,0,0])
at2 = Atom('O', [1.43, 0,0])
CO = at1 + at2
CO
Pt111 = surflab.fccsurfaces('Pt', '111', (2,2,4), vac=15)
grp = carbonlab.graphene(0,0)
from oorinano.visualizer import vis
vis.show_xcrysden(Pt111)
vis.show_xcrysden(grp)
