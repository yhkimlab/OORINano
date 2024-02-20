#
# Module Thermodynamics
#
# Developer   : Min Jong Noh
# Last update : 2020/10/01
# E-mail      : starnmj@kaist.ac.kr
# updated by J. Park 2021/10/18: refactoring: change class to module
# updated by J. Park 2023/01/18: modify

import operator
from ...atoms import *
from .calgibbs import pH_free_energy

def common_plot():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.use('Agg')
    plt.figure(figsize=(8.5, 6.5))
    plt.rcParams['axes.linewidth'] = 2.0

    return plt

def common_format(plt):
    ### common font for Gibbs plot
    plt.tick_params(axis="y", direction="in", length=10, width=2)
    plt.tick_params(axis="x", direction="in", length=10, width=2)
    plt.xlabel('Reaction Coordinate', fontsize=25)
    plt.ylabel('Gibbs Energy [eV]', fontsize=25)
    
    return 0


def plot_HER(Gibbs_H, legend=None, ymin=-1, ymax=1, label="HER", y_ticks=[-1.0, -0.5, 0, 0.5, 1.0], dpi=600, pH=0, Temp=298.15):
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

    plt = common_plot()

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
    
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={"size":15}, frameon=False)
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')

def plot_ORR_4e_wU(Gibbs_E, U=None, legend=None, ymin=None, ymax=None, label="ORR_4e", y_ticks=[-2, 0, 2, 4, 6], dpi=600, G_cal=0, pH=0, Temp=298.15):
    '''
    Plot Gibbs energy with U (cell potential)
        pH is included in Gibbs_E
        U is calculated in routine for U_oc (open-circuit, onset for overpotential) using U_0 and U_eq
    '''
    plt = common_plot()
    #global R, T
    X        = [-1, 0, 1, 2, 3, 4, 5, 6]
    # ymin=-2
    # ymax=6
    xmin=-0.5; xmax=5.5
    
    G_0     = Gibbs_E
    ### --> find onset potential for input
    pot_onset       = - max(list(map(operator.sub, Gibbs_E[1:],Gibbs_E[:-1]))) 
    print(f"onset potential {pot_onset}")

    # Gibbs free energy with applied potential
    G_pH = pH_free_energy(pH=pH, Temp=Temp)
    U_eq = 1.23 - G_pH
    eta=U_eq-pot_onset
    print(f"Equilibrium pot at pH{pH:2d} = {U_eq:.2f} with onset potential {pot_onset:.2f}, overpotential {eta:.2f}")

    applied_U = [0.00, pot_onset, U_eq]
    
    if legend is None:
        legend = [r'U$_0$=0.00 V', r'U$_{oc}$=%4.2f V' % pot_onset, r'U$_{eq}$=%4.2f V' % U_eq] # Voc == Vonset
    else:
        legend = legend
    label = label + f"_Voc{pot_onset:.1f}"
    
    elec_conv = -1
    X         = [-1, 0, 1, 2, 3, 4, 5, 6]
    G_eU=[]
    ### ref -> G[5], Note that pH is included in Gibbs in calc_gibbs_ORR_4e()
    for i in range(len(applied_U)):
        G_U   = [G_0[0]-G_0[5]+4*elec_conv*applied_U[i], # G_System+O2   + 4*e-, N.B pH is included in G_0
                 G_0[1]-G_0[5]+4*elec_conv*applied_U[i], # G_System-O2*  + 4*e-
                 G_0[2]-G_0[5]+3*elec_conv*applied_U[i], # G_System-OOH* + 3*e-
                 G_0[3]-G_0[5]+2*elec_conv*applied_U[i], # G_System+O    + 2*e-
                 G_0[4]-G_0[5]+1*elec_conv*applied_U[i], # G_System+OH   + 1*e-
                 G_0[5]-G_0[5]+0*elec_conv*applied_U[i], # G_System      + 0*e-
                ]
        G_U.insert(0, G_U[0])       # dup. for plot
        G_U.append(G_U[-1])         # dup. for plot
        G_eU.append(G_U)
        
    if G_cal:
        return G_U     
    
    #print(f"len(X) {len(X)}, len(Y) {len(G_eU[0])}")
    colors=['red', 'blue', 'black']
    for j in range(len(G_eU)):
        plt.scatter(X, G_eU[j], color=colors[j], marker='_', s=6000, linewidth=4, label='%s' % legend[j])
        for i in range(len(X)-1):
            plt.vlines((X[i+1]+X[i])/2, G_eU[j][i], G_eU[j][i+1], linestyle=':', color=colors[j], linewidth=2)
    
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4, 5], ['O$_{2}$', 'O$_{2}^{*}$', 'OOH$^{*}$', 'O$^{*}$', 'OH$^{*}$', 'H$_{2}$O'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    common_format(plt)
    #plt.legend(loc=2, prop={"size":10})
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), prop={"size":25}, frameon=False)
    
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')


def plot_OER_4e_wU(Gibbs_E, U=None, legend=None, ymin=None, ymax=None, label="OER_4e", y_ticks=[-2, 0, 2, 4, 6], dpi=600, G_cal=0, pH=0, Temp=298.15):
    '''
    Plot Gibbs energy with U (cell potential)
        pH is included in Gibbs_E
        U is calculated in routine for U_oc (open-circuit, onset for overpotential) using U_0 and U_eq
    '''
    plt = common_plot()
    #common_font(plt)
    X    = [-1, 0, 1, 2, 3, 4, 5]

    G_0  = Gibbs_E 
    
    pot_onset   = max(list(map(operator.sub, Gibbs_E[1:],Gibbs_E[:-1]))) 
    #print(f"onset potential {pot_onset:2.f}")

    # Gibbs free energy with applied potential
    G_pH = pH_free_energy(pH=pH, Temp=Temp)
    U_eq = 1.23 - G_pH
    eta = pot_onset - U_eq
    print(f"Equilibrium pot at pH{pH:2d} = {U_eq:.2f} with onset potential {pot_onset:.2f}, overpotential {eta:.2f}")


    applied_U = [0.00, U_eq, pot_onset]

    if legend is None:
        legend = [r'U$_0$=0.00 V', r'U$_{eq}$=%4.2f V' % U_eq, r'U$_{oc}$=%4.2f V' % pot_onset] # Voc == Vonset
    else:
        legend = legend
    label = label + f"_eta{eta:.1f}"

    elec_conv = -1
    
    G_eU=[]
    
    for i in range(len(applied_U)): # G_eU_max
        G_U  = [ G_0[0]-G_0[0]+0*elec_conv*applied_U[i], # G_System+2H2O + 0*e- N.B. pH is included in G_0
                 G_0[1]-G_0[0]+1*elec_conv*applied_U[i], # G_System-OH*  + 1*e-
                 G_0[2]-G_0[0]+2*elec_conv*applied_U[i], # G_System-O*   + 2*e-
                 G_0[3]-G_0[0]+3*elec_conv*applied_U[i], # G_System+OOH* + 3*e-
                 G_0[4]-G_0[0]+4*elec_conv*applied_U[i], # G_System+O2   + 4*e-, No O2* state
                 ]
        G_U.insert(0, G_U[0])       # dup. for plot
        G_U.append(G_U[-1])         # dup. for plot
        G_eU.append(G_U)

    if G_cal:
        return G_eU     
     
    colors=['red', 'black', 'blue']
    for j in range(len(G_eU)):
        plt.scatter(X, G_eU[j], color=colors[j], marker='_', s=9000, linewidth=4, label='%s' % legend[j])
        for i in range(len(X)-1):
            plt.vlines((X[i+1]+X[i])/2, G_eU[j][i], G_eU[j][i+1], linestyle=':', color=colors[j], linewidth=2)
                                                                                                                                 
    plt.axis([-0.50, 4.50, ymin, ymax])
    plt.xticks([0, 1, 2, 3, 4], ['H$_{2}$O', 'OH$^{*}$', 'O$^{*}$', 'OOH$^{*}$', 'O$_{2}$'], fontsize=20)
    plt.yticks(y_ticks, fontsize=20)
    common_format(plt)
    plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.5), prop={"size":25}, frameon=False)
    
    plt.savefig('%s.png' % label, format='png', dpi=dpi, bbox_inches = 'tight')
