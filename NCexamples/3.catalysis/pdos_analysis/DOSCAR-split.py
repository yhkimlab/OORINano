import NanoCore as nc

at = nc.io.read_poscar('POSCAR')
n_at = len(at)
nc.vasp.pdos_split_sum(sum_list=[i for i in range(n_at)])

