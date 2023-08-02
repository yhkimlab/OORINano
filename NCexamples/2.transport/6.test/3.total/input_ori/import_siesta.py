from nanocore.simulator import siesta as s2
from nanocore import *

model = s2.read_struct("STRUCT.fdf")
sim = s2.Siesta(model)
sim.read_fdf("BASIS.fdf")
sim.read_fdf("RUN.fdf")
sim.read_fdf("KPT.fdf")