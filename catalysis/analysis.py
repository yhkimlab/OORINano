#
# Module Thermodynamics
#
# Developer   : Min Jong Noh
# Last update : 2020/10/01
# E-mail      : starnmj@kaist.ac.kr
# updated by Joonho Park 2021/10/18: refactoring: change class to module

import os, sys, glob, math
from ..atoms import *
'''
def
    free_energies
    gibbs_HER
    gibbs_ORR_4e_acid
    gibbs_ORR_4e_alkaline
    gibbs_OER_4e_acid
'''
def free_energies(temp=298.15, pH=0, p=0.035, sol=0):

    E_H2O    = -14.236         # Total energy for H2O(g) in VASP
    if sol == 1:
        E_H2O = E_H2O - 0.310
    else:
        pass
    ZPE_H2O  =  0.560          # Norskov parameter
    S_H2O    =  0.00223333333  # Norskov parameter 0.67 eV at 300 K
    G_H2O_g = E_H2O + ZPE_H2O - temp*S_H2O
    
    E_H2   =  -6.760           # Total energy for H2(g) in VASP
    if sol == 1:
        E_H2 = E_H2 + 0.056
    else:
        pass
    ZPE_H2   = 0.270           # Norskov parameter
    TS_H2   =  0.001366666666  # Norskov parameter 0.41 eV at 300 K
    G_H2_g  = E_H2  + ZPE_H2  - temp*TS_H2
    
    R      = 0.0000861733254056734 # eV / K
    kB     = 0.0000861733254056734 # eV / K
    p0     = 1
    E_pH   = kB * temp * np.log(10) * pH
    
    # H2O
    G_H2O_aq = G_H2O_g + R*temp*np.log(p/p0) # only available for around p=0.035, T=298.15K

    # O2
    O2_g     = 2 * G_H2O_aq - 2 * G_H2_g + 4.92 # O2(g) from O2 + 2H2 -> 2H2O
    
    # H+
    H_ion    = 0.5*G_H2_g - E_pH
    
    # OH-
    OH_ion   = G_H2O_aq - H_ion
    
    return G_H2O_aq, O2_g, H_ion, OH_ion

def gibbs_HER(Sys, SysH, ZPE=None, TS=None):
    E_H2   =  -6.760;  ZPE_H2   = 0.270; TS_H2   = 0.410 
    
    n_component = len(Sys)

    if ZPE is None or TS is None:
        ZPE = []
        for i in range(n_component):
            ZPE_default = 0.04 + 0.5*ZPE_H2  # Cu111 + H (Norskov) J. Electro. Soc. 152(3), J23 (2005)
            ZPE.append(ZPE_default)
    else:
        ZPE = ZPE

    if TS is None:
        TS = []
        for i in range(n_component):
            TS_default  = 0.20 + 0.5*TS_H2   # Cu111 + H (Norskov)
            TS.append(TS_default)
    else:
        TS = TS

    if len(SysH) + len(ZPE) + len(TS) - 3*n_component == 0:
        pass
    else:
        print("The input components do not match")
        print("Check the input list")
        print("E_System:", len(Sys))
        print("E_System+H:", len(SysH))
        print("E_ZPE:", len(ZPE))
        print("E_TS:", len(TS))
   
    Gibbs_H = []

    for i in range(n_component):
        delta_E    = SysH[i] - (Sys[i] + 0.5*E_H2)
        delta_ZPE  = ZPE[i] - 0.5*ZPE_H2
        delta_TS   = TS[i]  - 0.5*TS_H2
        delta_G    = delta_E + delta_ZPE + delta_TS
        Gibbs_H.append(delta_G)

    return Gibbs_H

def gibbs_ORR_4e_acid(TE, ZPE=None, TS=None, temp=298.15, pH=0, p=0.035, sol=0):
    
    n_component = len(TE)  
    # TE (DFT total energy) must contain a series of energies
    # TE = [Sys, Sys+O2, Sys+OOH, Sys+O, Sys+OH]
                                                                                                                           
    if ZPE is None:
        ZPE = [0.000, 0.169, 0.511, 0.117, 0.430]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        ZPE = ZPE
    
    if TS is None:
        TS = [0.000, 0.097, 0.087, 0.020, 0.030]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        TS = TS

    if len(TE) + len(ZPE) + len(TS) - 3*n_component == 0:
        pass
    else:
        print("The input components do not match")
        print("Check the input list")
        print("TE:", len(TE))
        print("ZPE:", len(ZPE))
        print("TS:", len(TS))
     
    G_H2O_aq, O2_g, H_ion, OH_ion = free_energies(temp=temp, pH=pH, p=p, sol=sol)
    G_Sys     = TE[0] + ZPE[0] - TS[0] + 4 * H_ion + 0 * G_H2O_aq + 1 * O2_g
    G_SysO2   = TE[1] + ZPE[1] - TS[1] + 4 * H_ion + 0 * G_H2O_aq
    G_SysOOH  = TE[2] + ZPE[2] - TS[2] + 3 * H_ion + 0 * G_H2O_aq
    G_SysO    = TE[3] + ZPE[3] - TS[3] + 2 * H_ion + 1 * G_H2O_aq
    G_SysOH   = TE[4] + ZPE[4] - TS[4] + 1 * H_ion + 1 * G_H2O_aq                         
    G_Sys_end = TE[0] + ZPE[0] - TS[0] + 0 * H_ion + 2 * G_H2O_aq 
    Gibbs_E   = [G_Sys, G_SysO2, G_SysOOH, G_SysO, G_SysOH, G_Sys_end]

    return Gibbs_E

def gibbs_ORR_4e_alkaline(TE, ZPE=None, TS=None, temp=298.15, pH=14, p=0.035, sol=0):
    
    n_component = len(TE)  
    # TE (DFT total energy) must contain a series of energies
    # TE = [Sys, Sys+O2, Sys+OOH, Sys+O, Sys+OH]
                                                                                                                           
    if ZPE is None:
        ZPE = [0.000, 0.169, 0.511, 0.117, 0.430]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        ZPE = ZPE
    
    if TS is None:
        TS = [0.000, 0.097, 0.087, 0.020, 0.030]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        TS=TS
                                                                                                                   
    if len(TE) + len(ZPE) + len(TS) - 3*n_component == 0:
        pass
    else:
        print("The input components do not match")
        print("Check the input list")
        print("TE:", len(TE))
        print("ZPE:", len(ZPE))
        print("TS:", len(TS))
    
    G_H2O_aq, O2_g, H_ion, OH_ion = free_energies(temp=temp, pH=pH, p=p, sol=sol)
    G_Sys     = TE[0] + ZPE[0] - TS[0] + 0 * OH_ion + 2 * G_H2O_aq + 1 * O2_g
    G_SysO2   = TE[1] + ZPE[1] - TS[1] + 0 * OH_ion + 2 * G_H2O_aq
    G_SysOOH  = TE[2] + ZPE[2] - TS[2] + 1 * OH_ion + 1 * G_H2O_aq
    G_SysO    = TE[3] + ZPE[3] - TS[3] + 2 * OH_ion + 1 * G_H2O_aq
    G_SysOH   = TE[4] + ZPE[4] - TS[4] + 3 * OH_ion + 0 * G_H2O_aq                         
    G_Sys_end = TE[0] + ZPE[0] - TS[0] + 4 * OH_ion + 0 * G_H2O_aq 
    Gibbs_E   = [G_Sys, G_SysO2, G_SysOOH, G_SysO, G_SysOH, G_Sys_end]
    
    return Gibbs_E

def gibbs_OER_4e_acid(TE, ZPE=None, TS=None, temp=298.15, pH=0, p=0.035):
    
    n_component = len(TE)  
    # TE (DFT total energy) must contain a series of energies
    # TE = [Sys, Sys+O2, Sys+OOH, Sys+O, Sys+OH]
                                                                                                                           
    if ZPE is None:
        ZPE = [0.000, 0.430, 0.117, 0.511]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        ZPE = ZPE
    
    if TS is None:
        TS = [0.000, 0.030, 0.020, 0.087]
        # reference values from Nat.commun. 8, 15938 (2017) STable 7
    else:
        TS = TS
                                                                                                                            
    if len(TE) + len(ZPE) + len(TS) - 3*n_component == 0:
        pass
    else:
        print("The input components do not match")
        print("Check the input list")
        print("TE:", len(TE))
        print("ZPE:", len(ZPE))
        print("TS:", len(TS))
     
    G_H2O_aq, O2_g, H_ion, OH_ion = free_energies(temp=temp, pH=pH, p=p)
    G_Sys     = TE[0] + ZPE[0] - TS[0] + 0 * H_ion + 2 * G_H2O_aq 
    G_SysOH   = TE[1] + ZPE[1] - TS[1] + 1 * H_ion + 1 * G_H2O_aq
    G_SysO    = TE[2] + ZPE[2] - TS[2] + 2 * H_ion + 1 * G_H2O_aq
    G_SysOOH  = TE[3] + ZPE[3] - TS[3] + 3 * H_ion + 0 * G_H2O_aq
    G_Sys_end = TE[0] + ZPE[0] - TS[0] + 4 * H_ion + 0 * G_H2O_aq + 1 * O2_g 
    Gibbs_E   = [G_Sys, G_SysOH, G_SysO, G_SysOOH, G_Sys_end]
                                                                                                                            
    return Gibbs_E

