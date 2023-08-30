#
# Module Thermodynamics
#
# Developer   : Min Jong Noh
# Last update : 2020/10/01
# E-mail      : starnmj@kaist.ac.kr
# updated by J. Park 2021/10/18: refactoring: change class to module
# updated by J. Park 2023/01/18: modify

import os, sys, glob, math, operator
from ...atoms import *
from .calgibbs import pH_free_energy


def plot_HER(Gibbs_H, legend=None, ymin=-1, ymax=1, label="HER", y_ticks=[-1.0, -0.5, 0, 0.5, 1.0], dpi=600):
    """
    plot HER profile with python-matplotlib
    this module is independent, so if you can use only hand-shaved values
    
    Parameters
    ----------
    Gibbs_H: list
       Gibbs free energy which to be directly plotted.
       Ex. [-0.32, -0.11, 0.12]
       tips: module run_HER returns Gibbs_H value, therefore, append the value on the list.
    
    legend:
       the legend when matplotlib will visualize.
       this is only available for number of dataset is matched between Gibbs_H and legend.
       otherwise, the legend is automatically assigend just [1, 2, 3 ...]
    """
    e_Gibbs_H = Gibbs_H

    X = [0, 1, 2]
    Z = []

    if legend is None:
        for i in range(len(e_Gibbs_H)):
            Z.append(i+1)
    elif len(Gibbs_H) == len(legend):
        Z = legend
    else:
        raise ValueError("The number of data between Gibbs energies and their legends is not matched")

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, 8.5))
    plt.rcParams['axes.linewidth'] = 2.0
    X_domain    = [-0.5,  0.0,  2.5]
    upper_limit = [ 0.2,  0.2,  0.2]
    lower_limit = [-0.2, -0.2, -0.2]

    hori1 = plt.plot(X_domain, upper_limit, linestyle='dotted', linewidth='0.5', color='k')
    hori2 = plt.plot(X_domain, lower_limit, linestyle='dotted', linewidth='0.5', color='k')
    plt.fill_between(X_domain, upper_limit, 0, color = 'whitesmoke')
    plt.fill_between(X_domain, lower_limit, 0, color = 'whitesmoke')

    for i in range(len(e_Gibbs_H)):
        Y = [0, e_Gibbs_H[i], 0]
        plt.plot(X, Y, marker='_', ms='50', mew='5', linestyle='dashed', linewidth='1', label='%s' % str(Z[i]))
    plt.axis([-0.4, 2.4, ymin, ymax])
    plt.xticks([0, 1, 2], [r' H$^{+}$+e$^{-}$ ', r' H$^{*}$ ', r' 1/2H$_{2}$ '], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Free Energy [eV]', fontsize=25)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={"size":15}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')

def plot_ORR_4e(Gibbs_E, U=None, legend=None, ymin=None, ymax=None, label="ORR_4e", y_ticks=[-2, 0, 2, 4, 6], dpi=600, G_cal=0, pH=0, T=298.15):
    #global R, T
    X        = [-1, 0, 1, 2, 3, 4, 5, 6]
    # ymin=-2
    # ymax=6
    xmin=-0.5; xmax=5.5
    
    #G_eq  = Gibbs_E 
    G_0     = Gibbs_E
    ### --> find onset potential for input
    pot_onset       = - max(list(map(operator.sub, Gibbs_E[1:],Gibbs_E[:-1]))) 
    print(pot_onset)

    # Gibbs free energy with applied potential
    G_pH = pH_free_energy(pH=pH, T=T)
    U_eq = 1.23 - G_pH
    print(f"Equilibrium pot at pH{pH:2d} = {U_eq} with onset potential {pot_onset:.2f}, overpotential {U_eq-pot_onset:.2f}")

    applied_U = [0.00, pot_onset, U_eq]
    
    if legend is None:
        legend = [r'U$_0$=0.00 V', r'U$_{onset}$=%4.2f V' % pot_onset, r'U$_{eq}$=%4.2f V' % U_eq]
    else:
        legend = legend
    label = label + f"_eta{pot_onset:.1f}"
    
    elec_conv = -1
    X         = [-1, 0, 1, 2, 3, 4, 5, 6]
    G_eU=[]
    for i in range(len(applied_U)):
        G_U   = [G_0[0]-G_0[5]+4*elec_conv*applied_U[i], # dummy for first step plot
                 G_0[0]-G_0[5]+4*elec_conv*applied_U[i], # G_System+O2   + 4*e-
                 G_0[1]-G_0[5]+4*elec_conv*applied_U[i], # G_System-O2*  + 4*e-
                 G_0[2]-G_0[5]+3*elec_conv*applied_U[i], # G_System-OOH* + 3*e-
                 G_0[3]-G_0[5]+2*elec_conv*applied_U[i], # G_System+O    + 2*e-
                 G_0[4]-G_0[5]+1*elec_conv*applied_U[i], # G_System+OH   + 1*e-
                 G_0[5]-G_0[5]+0*elec_conv*applied_U[i], # G_System      + 0*e-
                 G_0[5]-G_0[5]+0*elec_conv*applied_U[i], # dummy for last step plot
                ]
        G_eU.append(G_U)

    if G_cal:
        return G_U     
    else:
        pass
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8.5, 6.5))
    plt.rcParams['axes.linewidth'] = 2.0
    
    X_dom       = [-2,  0, 4, 10]
    Y_hori      = [ 0,  0, 0,  0]
    plt.scatter(X, G_eU[0], color='red', marker='_', s=6000, linewidth=4, label='%s' % legend[0])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[0][i], G_eU[0][i+1], linestyle=':', color='red', linewidth=2)
    plt.scatter(X, G_eU[1], color='blue', marker='_', s=6000, linewidth=4, label='%s' % legend[1])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[1][i], G_eU[1][i+1], linestyle=':', color='blue', linewidth=2)
    plt.scatter(X, G_eU[2], color='black', marker='_', s=6000, linewidth=4, label='%s' % legend[2])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[2][i], G_eU[2][i+1], linestyle=':', color='black', linewidth=2)           
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4, 5], ['O$_{2}$', 'O$_{2}^{*}$', 'OOH$^{*}$', 'O$^{*}$', 'OH$^{*}$', 'H$_{2}$O'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Free Energy [eV]', fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), prop={"size":25}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')



def plot_ORR_4e_acid(Gibbs_E, U=None, legend=None, ymin=None, ymax=None, label="ORR_4e", y_ticks=[-2, 0, 2, 4, 6], dpi=600, G_cal=0):
    
    X        = [-1, 0, 1, 2, 3, 4, 5, 6]
    # ymin=-2
    # ymax=6
    xmin=-0.5; xmax=5.5
    G_eq  = Gibbs_E 
    
    # Gibbs free energy with applied potential
    
    applied_U = [1.23]
    
    elec_conv = -1
    X         = [-1, 0, 1, 2, 3, 4, 5, 6]
    G_eU_max  = [G_eq[0]-G_eq[5]+4*elec_conv*applied_U[0], # dummy for first step plot
                 G_eq[0]-G_eq[5]+4*elec_conv*applied_U[0], # G_System+O2   + 4*e-
                 G_eq[1]-G_eq[5]+4*elec_conv*applied_U[0], # G_System-O2*  + 4*e-
                 G_eq[2]-G_eq[5]+3*elec_conv*applied_U[0], # G_System-OOH* + 3*e-
                 G_eq[3]-G_eq[5]+2*elec_conv*applied_U[0], # G_System+O    + 2*e-
                 G_eq[4]-G_eq[5]+1*elec_conv*applied_U[0], # G_System+OH   + 1*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[0], # G_System      + 0*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[0], # dummy for last step plot
                ]

    G_rel     = []
  
    for i in range(len(G_eU_max)-1):      
        rel_G = G_eU_max[i+1] - G_eU_max[i]   
        G_rel.append(rel_G)           
                                    
    if U is None:                    
        U_theory = 1.23 - max(G_rel) 
        if U_theory < 0:
            U = 0
        else:
            U = U_theory                 
    else:                            
        U = U                        
    
    if legend is None:
        legend = ['U=1.23V', 'U=%3.2fV' % U, 'U=0.00V']
    else:
        legend = legend
    
    applied_U = [1.23, U, 0]
    
    G_eU      = [G_eq[0]-G_eq[5]+4*elec_conv*applied_U[1], # dummy for first step plot
                 G_eq[0]-G_eq[5]+4*elec_conv*applied_U[1], # G_System+O2    + 4*e-
                 G_eq[1]-G_eq[5]+4*elec_conv*applied_U[1], # G_System-O2*   + 4*e-
                 G_eq[2]-G_eq[5]+3*elec_conv*applied_U[1], # G_System-OOH*  + 3*e-
                 G_eq[3]-G_eq[5]+2*elec_conv*applied_U[1], # G_System+O     + 2*e-
                 G_eq[4]-G_eq[5]+1*elec_conv*applied_U[1], # G_System+OH    + 1*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[1], # G_System       + 0*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[1], # dummy for last step plot
                ]
    
    G_eU_min  = [G_eq[0]-G_eq[5]+4*elec_conv*applied_U[2], # dummy for first step plot
                 G_eq[0]-G_eq[5]+4*elec_conv*applied_U[2], # G_System+O2    + 4*e-
                 G_eq[1]-G_eq[5]+4*elec_conv*applied_U[2], # G_System-O2*   + 4*e-
                 G_eq[2]-G_eq[5]+3*elec_conv*applied_U[2], # G_System-OOH*  + 3*e-
                 G_eq[3]-G_eq[5]+2*elec_conv*applied_U[2], # G_System+O     + 2*e-
                 G_eq[4]-G_eq[5]+1*elec_conv*applied_U[2], # G_System+OH    + 1*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[2], # G_System       + 0*e-
                 G_eq[5]-G_eq[5]+0*elec_conv*applied_U[2], # dummy for last step plot
               ]
    
    if G_cal:
        return G_eU     
    else:
        pass
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8.5, 6.5))
    plt.rcParams['axes.linewidth'] = 2.0
    
    X_dom       = [-2,  0, 4, 10]
    Y_hori      = [ 0,  0, 0,  0]
    plt.scatter(X, G_eU_max, color='black', marker='_', s=6000, linewidth=4, label='%s' % legend[0])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_max[i], G_eU_max[i+1], linestyle=':', color='black', linewidth=2)
    plt.scatter(X, G_eU, color='red', marker='_', s=6000, linewidth=4, label='%s' % legend[1])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[i], G_eU[i+1], linestyle=':', color='red', linewidth=2)
    plt.scatter(X, G_eU_min, color='blue', marker='_', s=6000, linewidth=4, label='%s' % legend[2])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_min[i], G_eU_min[i+1], linestyle=':', color='blue', linewidth=2)           
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4, 5], ['O$_{2}$', 'O$_{2}^{*}$', 'OOH$^{*}$', 'O$^{*}$', 'OH$^{*}$', 'H$_{2}$O'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Free Energy [eV]', fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), prop={"size":25}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')

def plot_ORR_4e_alkaline(Gibbs_E, U=None, legend=None, ymin=-6, ymax=2, label="ORR_4e_alkaline", y_ticks=[-6,-4,-2,0,2],dpi=600, G_cal=0):
    
    X        = [-1, 0, 1, 2, 3, 4, 5, 6]
    
    G_eq  = Gibbs_E 
    
    # Gibbs free energy with applied potential
    
    applied_U = [0.40, U, -0.83]
    elec_conv = -1
    X         = [-1, 0, 1, 2, 3, 4, 5, 6]
    G_eU_max_ref = G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[0] 
    G_eU_max  = [G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[0]-G_eU_max_ref, # dummy for first step plot
                 G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+O2   + 4*e-
                 G_eq[1]-G_eq[5]+(4)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System-O2*  + 4*e-
                 G_eq[2]-G_eq[5]+(3)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System-OOH* + 3*e-
                 G_eq[3]-G_eq[5]+(2)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+O    + 2*e-
                 G_eq[4]-G_eq[5]+(1)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+OH   + 1*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[0]-G_eU_max_ref, # G_System      + 0*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[0]-G_eU_max_ref, # dummy for last step plot
                ]
                                                                                                                            
    G_rel     = []
                                                        
    for i in range(len(G_eU_max)-1):      
        rel_G = G_eU_max[i+1] - G_eU_max[i]   
        G_rel.append(rel_G)           
                                    
    if U is None:                    
        U_theory = 0.40 - max(G_rel) 
        if U_theory < -0.83:
            U = -0.83
        else:
            U = U_theory                 
    else:                            
        U = U                        
    
    if legend is None:
        legend = ['U=0.40V', 'U=%3.2fV' % U, 'U=-0.83V']
    else:
        legend = legend
    
    applied_U = [0.40, U, -0.83]
    G_eU_ref  = G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[1]
    G_eU      = [G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[1]-G_eU_ref, # dummy for first step plot
                 G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[1]-G_eU_ref, # G_System+O2    + 4*e-
                 G_eq[1]-G_eq[5]+(4)*elec_conv*applied_U[1]-G_eU_ref, # G_System-O2*   + 4*e-
                 G_eq[2]-G_eq[5]+(3)*elec_conv*applied_U[1]-G_eU_ref, # G_System-OOH*  + 3*e-
                 G_eq[3]-G_eq[5]+(2)*elec_conv*applied_U[1]-G_eU_ref, # G_System+O     + 2*e-
                 G_eq[4]-G_eq[5]+(1)*elec_conv*applied_U[1]-G_eU_ref, # G_System+OH    + 1*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[1]-G_eU_ref, # G_System       + 0*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[1]-G_eU_ref, # dummy for last step plot
                ]
    G_eU_min_ref = G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[2]
    G_eU_min  = [G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[2]-G_eU_min_ref, # dummy for first step plot
                 G_eq[0]-G_eq[5]+(4)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+O2    + 4*e-
                 G_eq[1]-G_eq[5]+(4)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System-O2*   + 4*e-
                 G_eq[2]-G_eq[5]+(3)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System-OOH*  + 3*e-
                 G_eq[3]-G_eq[5]+(2)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+O     + 2*e-
                 G_eq[4]-G_eq[5]+(1)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+OH    + 1*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[2]-G_eU_min_ref, # G_System       + 0*e-
                 G_eq[5]-G_eq[5]+(0)*elec_conv*applied_U[2]-G_eU_min_ref, # dummy for last step plot
               ]
    
    if G_cal:
        return G_eU    
    else:
        pass

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8.5, 6.5))
    plt.rcParams['axes.linewidth'] = 2.0
    
    X_dom       = [-2,  0, 4, 10]
    Y_hori      = [ 0,  0, 0,  0]
    plt.scatter(X, G_eU_max, color='black', marker='_', s=6000, linewidth=4, label='%s' % legend[0])                 
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_max[i], G_eU_max[i+1], linestyle=':', color='black', linewidth=2)
    plt.scatter(X, G_eU, color='red', marker='_', s=6000, linewidth=4, label='%s' % legend[1])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[i], G_eU[i+1], linestyle=':', color='red', linewidth=2)
    plt.scatter(X, G_eU_min, color='blue', marker='_', s=6000, linewidth=4, label='%s' % legend[2])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_min[i], G_eU_min[i+1], linestyle=':', color='blue', linewidth=2)           
                                                                                                                            
    plt.axis([-0.5, 5.5, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4, 5], ['O$_{2}$', 'O$_{2}^{*}$', 'OOH$^{*}$', 'O$^{*}$', 'OH$^{*}$', 'OH$^{-}$'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Free Energy [eV]', fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), prop={"size":25}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')

def plot_OER_4e_acid(Gibbs_E, U=None, legend=None, ymin=-6, ymax=6, label="OER_4e_acid", y_ticks=[-6, -4, -2, 0, 2, 4, 6], dpi=600, G_cal=0):
    
    X        = [-1, 0, 1, 2, 3, 4, 5]
    
    G_eq  = Gibbs_E 
    
    # Gibbs free energy with applied potential
    
    applied_U = [1.23]
    
    elec_conv = -1
    X         = [-1, 0, 1, 2, 3, 4, 5]
    G_eU_max_ref = G_eq[0]-G_eq[4]+0*elec_conv*applied_U[0]
    G_eU_max  = [G_eq[0]-G_eq[4]+0*elec_conv*applied_U[0]-G_eU_max_ref, # dummy for first step plot
                 G_eq[0]-G_eq[4]+0*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+2H2O + 0*e-
                 G_eq[1]-G_eq[4]+1*elec_conv*applied_U[0]-G_eU_max_ref, # G_System-OH*  + 1*e-
                 G_eq[2]-G_eq[4]+2*elec_conv*applied_U[0]-G_eU_max_ref, # G_System-O*   + 2*e-
                 G_eq[3]-G_eq[4]+3*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+OOH* + 3*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[0]-G_eU_max_ref, # G_System+O2   + 4*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[0]-G_eU_max_ref, # dummy for last step plot 
                ]
                                                                                                                                            
    G_rel     = []
  
    for i in range(len(G_eU_max)-1):      
        rel_G = G_eU_max[i+1] - G_eU_max[i]   
        G_rel.append(rel_G)           
                                    
    if U is None:                    
        U_theory = 1.23 + max(G_rel) 
        if U_theory > 2.46:
            U = 2.46
        else:
            U = U_theory                 
    else:                            
        U = U                        
    
    if legend is None:
        legend = ['U=1.23V', 'U=%3.2fV' % U, 'U=0.00V']
    else:
        legend = legend
    
    applied_U = [1.23, U, 0]
    G_eU_ref  =  G_eq[0]-G_eq[4]+0*elec_conv*applied_U[1],
    G_eU      = [G_eq[0]-G_eq[4]+0*elec_conv*applied_U[1]-G_eU_ref, # dummy for first step plot
                 G_eq[0]-G_eq[4]+0*elec_conv*applied_U[1]-G_eU_ref, # G_System+2H2O  + 0*e-
                 G_eq[1]-G_eq[4]+1*elec_conv*applied_U[1]-G_eU_ref, # G_System-OH*   + 1*e-
                 G_eq[2]-G_eq[4]+2*elec_conv*applied_U[1]-G_eU_ref, # G_System-O*    + 2*e-
                 G_eq[3]-G_eq[4]+3*elec_conv*applied_U[1]-G_eU_ref, # G_System+OOH*  + 3*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[1]-G_eU_ref, # G_System+O2    + 4*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[1]-G_eU_ref, # dummy for last step plot
                ]
    
    G_eU_min_ref = G_eq[0]-G_eq[4]+0*elec_conv*applied_U[2],
    G_eU_min  = [G_eq[0]-G_eq[4]+0*elec_conv*applied_U[2]-G_eU_min_ref, # dummy for first step plot
                 G_eq[0]-G_eq[4]+0*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+2H2O  + 0*e-
                 G_eq[1]-G_eq[4]+1*elec_conv*applied_U[2]-G_eU_min_ref, # G_System-OH*   + 1*e-
                 G_eq[2]-G_eq[4]+2*elec_conv*applied_U[2]-G_eU_min_ref, # G_System-O*    + 2*e-
                 G_eq[3]-G_eq[4]+3*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+OOH*  + 3*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[2]-G_eU_min_ref, # G_System+O2    + 4*e-
                 G_eq[4]-G_eq[4]+4*elec_conv*applied_U[2]-G_eU_min_ref, # dummy for last step plot
               ]
    
    if G_cal:
        return G_eU     
    else:
        pass
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8.5, 6.5))
    plt.rcParams['axes.linewidth'] = 2.0
    
    X_dom       = [-2,  0, 4, 10]
    Y_hori      = [ 0,  0, 0,  0]
    plt.scatter(X, G_eU_max, color='black', marker='_', s=8000, linewidth=4, label='%s' % legend[0])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_max[i], G_eU_max[i+1], linestyle=':', color='black', linewidth=2)
    plt.scatter(X, G_eU, color='red', marker='_', s=8000, linewidth=4, label='%s' % legend[1])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU[i], G_eU[i+1], linestyle=':', color='red', linewidth=2)
    plt.scatter(X, G_eU_min, color='blue', marker='_', s=8000, linewidth=4, label='%s' % legend[2])
    for i in range(len(X)-1):
        plt.vlines((X[i+1]+X[i])/2, G_eU_min[i], G_eU_min[i+1], linestyle=':', color='blue', linewidth=2)           
                                                                                                                                            
    plt.axis([-0.50, 4.50, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4], ['H$_{2}$O', 'OH$^{*}$', 'O$^{*}$', 'OOH$^{*}$', 'O$_{2}$'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Free Energy [eV]', fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), prop={"size":25}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')
