# coding: utf-8
from oorinano import surflab
from oorinano.simulator.catalysis.catmodels import Catmodels
from oorinano.visualizer import vis
Pt111 = surflab.fccsurfaces('Pt', '111', (3,3,3), vac=15)
Pt_O2, Pt_OOH, Pt_O, Pt_OH = Catmodels(Pt111).four_electron_intermediates_gen(mode='ORR')
vis.show_xcrysden(Pt_OOH)
