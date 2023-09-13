import nanocore.io as nc

atms = nc.read_struct('elecL.fdf')
atms.get_mirrored_structure()
nc.write_struct(atms)