import os, sys, shutil
import numpy as np
import yaml

import matplotlib.pyplot as plt
from nanocore.simulators import siesta as s2
from mpl_toolkits.axes_grid1 import make_axes_locatable
import importlib


def qtPlot(calc, plot_dir, cal_dir, inpf, outpf='result', energy=[-0.1, 0.1], grid=[50,50,200]):
    simmodule = importlib.import_module(calc.__class__.__module__)

    cwd = os.getcwd()
    print(f"{cwd} in qtPLot")
    
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)
        print(f"{plot_dir} was generated")

    shutil.copy(f'{cwd}/{cal_dir}/{inpf}', f'{cwd}/{plot_dir}/')

    os.chdir(plot_dir)
    cwdw = os.getcwd()
    
    if 'files' not in os.listdir():
        os.mkdir('files')
    
    with open(inpf) as f:
        option = yaml.safe_load(f)
    fcwd = os.path.join(cwdw, 'files')
    os.chdir(option['scatter'])
    files = [f for f in os.listdir() if os.path.splitext(f)[-1] == '.VH']
    assert len(files) == 1
    shutil.copy(files[0], fcwd)
    os.chdir(fcwd)

    pos, pot = simmodule.get_hartree_pot_z(os.path.splitext(files[0])[0])
    pos = np.asarray(pos); pos = pos * 0.5291
    fig = plt.figure(figsize=[10,8])
    ax1 = fig.add_subplot(2,2,1)
    ax1.set_title("Potential", fontsize = 15)
    ax1.plot(pos, pot,'k', linewidth=3)
    ax1.set_xlabel('Z [Ang]',fontsize = 12)
    ax1.set_ylabel('$V_{H} [eV]$',fontsize = 12)
    ax1.tick_params(axis='both' ,labelsize=12, length = 6, width = 3, direction='in')
    plt.setp(ax1.spines.values(), linewidth=4)
    [ax1.spines[i].set_linewidth(3) for i in ax1.spines.keys()]
    ax1.tick_params(axis='both',which='major', labelsize=12,direction='in', length = 12, width = 3)
    ax1.tick_params(axis='both',which='minor', labelsize=12,direction='in', length = 6, width = 3)
    ax1.ticklabel_format(useOffset=False)
    fig.tight_layout()

    os.chdir(option['tbtrans'])
    files = os.listdir()
    for f in files:
        if os.path.splitext(f)[-1] in ['.AVTRANS_Left-Right']:
            file = f
            shutil.copy(f, fcwd)
    os.chdir(fcwd)

    energy, trans = simmodule.get_transmission(file)
    # trans = np.log10(trans)

    ax2 = fig.add_subplot(2,2,2)
    ax2.plot(energy, trans,'k', linewidth=3)
    ax2.set_title("Transmission", fontsize = 15)
    ax2.set_yscale('log')
    ax2.set_xlabel('$E-E_F [eV]$',fontsize = 12)
    ax2.set_ylabel('T',fontsize = 12)
    # ax1.set_ylim(0.01, 1)
    ax2.tick_params(axis='both' ,labelsize=12, length = 6, width = 3, direction='in')
    plt.setp(ax2.spines.values(), linewidth=4)
    [ax2.spines[i].set_linewidth(3) for i in ax1.spines.keys()]
    ax2.tick_params(axis='both',which='major', labelsize=12,direction='in', length = 12, width = 3)
    ax2.tick_params(axis='both',which='minor', labelsize=12,direction='in', length = 6, width = 3)
    fig.tight_layout()

    z_coords, energy, dos = s2.get_tbtrans_pldos(**option)

    X, Y = np.meshgrid(z_coords, energy)
    ax3 = fig.add_subplot(2,2,3)
    ax3.set_title("PLDOS", fontsize = 15)
    ax3.tick_params(axis='y', direction='in', length=12, width = 3,  pad=6, labelsize=12, labelcolor='black', right=False)
    ax3.tick_params(axis='x', direction='in', length=12, width = 3,  pad=6, labelsize=12, labelcolor='black', top=False)
    ax3.set_xlabel('$z [\AA]$',fontsize = 12)
    ax3.set_ylabel('$E-E_F [eV]$',fontsize = 12)

    levels = np.linspace(1.01*dos.min(), 0.99*dos.max(), 500)
    cset = ax3.contourf(X, Y, dos, levels, cmap='jet')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cset, cax=cax, orientation='vertical')
    fig.tight_layout()

    simmodule.get_tbtrans_ldos(os.path.join(cwdw, f'{outpf}.cube'), energy, grid, **option)

    os.chdir(cwdw)
    plt.savefig(f'{outpf}.jpg')
    os.chdir(cwd)
