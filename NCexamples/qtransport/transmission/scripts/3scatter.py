from nanocore.simulator import siesta as s2
import shutil
import os, sys
import yaml
import argparse

def run_qt(model, fin, fout, nproc):
    cwd = os.getcwd()
    #print(f"in 3scatter cwd {cwd}")
    shutil.copy(f'../2.model/{model}', 'input/STRUCT.fdf')      # copy model from 2.model/
    with open(fin) as f:
        yoption = yaml.safe_load(f)

    os.chdir('input')
    sim = s2.Siesta()
    sim.read_all_fdf()
    os.chdir('..')

    voltage = yoption['Voltage'].split()[0]
    if voltage not in os.listdir():
        os.mkdir(voltage)
        os.chdir(voltage)
        cwdv = os.getcwd()

        shutil.copytree('../input', 'TSHS')
        os.chdir('TSHS')
        sim.set_mode('scatter')
        sim.set_option('TS.Elecs.Neglect.Principal', True)
        sim.run_qt(nproc, **yoption)
        yoption['scatter'] = os.getcwd()
        os.chdir(cwdv)

        shutil.copytree('../input', 'TBTrans')
        os.chdir('TBTrans')
        sim.set_mode('tbtrans')
        sim.read_all_fdf()
        sim.set_option('kgrid_Monkhorst_Pack', True, (12,1,1))
        sim.set_option('TS.Elecs.Neglect.Principal', True)
        sim.run_qt(nproc, **yoption)
        yoption['tbtrans'] = os.getcwd()
        os.chdir(cwd)
        with open(fout, 'w') as f:
            yaml.safe_dump(yoption, f)
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",    default="input.yaml")
    parser.add_argument("-o", "--output",   default="output.yaml")
    parser.add_argument("-m", "--model",    default="cnt_6.fdf")
    parser.add_argument("-n", "--nproc", type=int, default=1)
    args = parser.parse_args()
    
    run_qt(args.model, args.input, args.output, args.nproc)
    return 0

if __name__ == '__main__':
    main()
