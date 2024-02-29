# coding: utf-8
from oorinano import surflab
from oorinano.simulator.catalysis.catmodels import Catmodels
from oorinano.visualizer import vis
Pt100 = surflab.fccsurfaces('Pt', '100', (2,2,4), vac=15)
Pt100H=Catmodels(Pt100).HER_intermediate_gen()
vis.show_xcrysden(Pt100H)
