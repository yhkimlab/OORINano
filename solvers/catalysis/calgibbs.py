#
# Module Thermodynamics
#
# Developer   : Min Jong Noh
# Last update : 2020/10/01
# E-mail      : starnmj@kaist.ac.kr
# updated by Joonho Park 2021/10/18: refactoring: change class to module

import os, sys, glob, math
from ...atoms import *
from ...thermo import R, p0, Etot_H2O, zpe_H2O, S_H2O, Etot_H2, zpe_H2, S_H2

def pH_free_energy(pH=0, T=298.15):
     
    G_pH    = R * T * np.log(10) * pH
    return G_pH


def mol_free_energies(T=298.15, pH=0, p=0.035, sol=0):
    '''
    Absolute Gibbs energy calculation for molecule and water
    '''
    #global Etot_H2O, zpe_H2O, Etot_H2, zpe_H2, R, kB, p0
    
    #if sol == 1:
    #    Etot_H2O -= 0.310
    #    Etot_H2 += 0.056
    G_H2O_g = Etot_H2O + zpe_H2O - T*S_H2O
    
    G_H2_g  = Etot_H2  + zpe_H2  - T*S_H2           # Gas phase calculation
    
    G_pH    = pH_free_energy(pH=0, T=298.15)        # RT ln10 = 0.05916
    
    ### Gibbs correction for H2O(l) from H2O(g, DFT)
    G_H2O_l = G_H2O_g + R*T*np.log(p/p0) # only available for around p=0.035,(0.0313?) T=298.15K

    # O2
    G_O2_g    = 2 * G_H2O_l - 2 * G_H2_g + 4.92     # O2(g) from O2 + 2H2 -> 2H2O
    
    # H+
    G_H_ion   = 0.5*G_H2_g - G_pH                   ### delG = delG0 - 0.0592pH
    
    # OH-
    G_OH_ion  = G_H2O_l - G_H_ion
    #print(f"G_H2O(l) {G_H2O_l:10.3f}, G_O2(g) {G_O2_g:10.3f}, G_H+ {G_H_ion:10.3f}, G_OH- {G_OH_ion:10.3f}")
    return G_H2O_l, G_O2_g, G_H_ion, G_OH_ion

def calc_gibbs_HER(Sys, SysH, ZPE=None, TS=None, Temp=298.15):
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

def calc_gibbs_ORR_4e(totE, ZPE=None, TS=None, Temp=298.15, pH=0, p=0.035, sol=0):
    
    n_component = len(totE)  
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

    if len(totE) + len(ZPE) + len(TS) - 3*n_component == 0:
        pass
    else:
        print("The input components do not match")
        print("Check the input list")
        print("TE:", len(totE))
        print("ZPE:", len(ZPE))
        print("TS:", len(TS))
    ### Gibbs for adsorbates 
    print("get Gibbs enegy for molecules")
    G_H2O_l, G_O2_g, G_H_ion, G_OH_ion = mol_free_energies(T=Temp, pH=pH, p=p, sol=sol)
    ### Gibbs through ORR
    print("calculate gibbs")
    G_Sys     = totE[0] + ZPE[0] - TS[0] + 4 * G_H_ion + 0 * G_H2O_l + 1 * G_O2_g
    G_SysO2   = totE[1] + ZPE[1] - TS[1] + 4 * G_H_ion + 0 * G_H2O_l
    G_SysOOH  = totE[2] + ZPE[2] - TS[2] + 3 * G_H_ion + 0 * G_H2O_l
    G_SysO    = totE[3] + ZPE[3] - TS[3] + 2 * G_H_ion + 1 * G_H2O_l
    G_SysOH   = totE[4] + ZPE[4] - TS[4] + 1 * G_H_ion + 1 * G_H2O_l                         
    G_Sys_end = totE[0] + ZPE[0] - TS[0] + 0 * G_H_ion + 2 * G_H2O_l 
    Gibbs_E   = [G_Sys, G_SysO2, G_SysOOH, G_SysO, G_SysOH, G_Sys_end]

    return Gibbs_E

def calc_gibbs_ORR_4e_alkaline(TE, ZPE=None, TS=None, T=298.15, pH=14, p=0.035, sol=0):
    
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
    
    G_H2O_l, O2_g, H_ion, OH_ion, G_pH = mol_free_energies(T=T, pH=pH, p=p, sol=sol)
    G_Sys     = TE[0] + ZPE[0] - TS[0] + 0 * OH_ion + 2 * G_H2O_l + 1 * O2_g
    G_SysO2   = TE[1] + ZPE[1] - TS[1] + 0 * OH_ion + 2 * G_H2O_l
    G_SysOOH  = TE[2] + ZPE[2] - TS[2] + 1 * OH_ion + 1 * G_H2O_l
    G_SysO    = TE[3] + ZPE[3] - TS[3] + 2 * OH_ion + 1 * G_H2O_l
    G_SysOH   = TE[4] + ZPE[4] - TS[4] + 3 * OH_ion + 0 * G_H2O_l                         
    G_Sys_end = TE[0] + ZPE[0] - TS[0] + 4 * OH_ion + 0 * G_H2O_l 
    Gibbs_E   = [G_Sys, G_SysO2, G_SysOOH, G_SysO, G_SysOH, G_Sys_end]
    
    return Gibbs_E, G_pH

def calc_gibbs_OER_4e_acid(TE, ZPE=None, TS=None, T=298.15, pH=0, p=0.035):
    
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
     
    G_H2O_l, O2_g, H_ion, OH_ion, G_pH = mol_free_energies(T=T, pH=pH, p=p)
    G_Sys     = TE[0] + ZPE[0] - TS[0] + 0 * H_ion + 2 * G_H2O_l 
    G_SysOH   = TE[1] + ZPE[1] - TS[1] + 1 * H_ion + 1 * G_H2O_l
    G_SysO    = TE[2] + ZPE[2] - TS[2] + 2 * H_ion + 1 * G_H2O_l
    G_SysOOH  = TE[3] + ZPE[3] - TS[3] + 3 * H_ion + 0 * G_H2O_l
    G_Sys_end = TE[0] + ZPE[0] - TS[0] + 4 * H_ion + 0 * G_H2O_l + 1 * O2_g 
    Gibbs_E   = [G_Sys, G_SysOH, G_SysO, G_SysOOH, G_Sys_end]
                                                                                                                            
    return Gibbs_E, G_pH

