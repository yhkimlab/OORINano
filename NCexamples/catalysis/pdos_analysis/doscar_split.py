import nanocore as nc

at = nc.vasp.read_poscar('POSCAR')
n_at = len(at)
nc.vasp.pdos_split_sum(sum_list=[i for i in range(n_at)])

