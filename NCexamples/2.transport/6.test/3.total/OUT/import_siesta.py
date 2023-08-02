from nanocore.simulator import siesta as s2
import yaml

with open("input.yaml") as f:
    option = yaml.safe_load(f)
model = s2.read_struct("STRUCT.fdf")
sim = s2.Siesta(model)
sim.read_fdf("BASIS.fdf")
sim.read_fdf("RUN.fdf")
sim.read_fdf("KPT.fdf")
sim.read_fdf("TS.fdf")
sim.run("scatter", 12, **option)
