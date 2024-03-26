
import os
#, shutil
#import argparse
#import sys
#import time

import argparse

# QFL-----------------------------------------------------------------------------
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection, PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
import glob, math
from matplotlib.cbook import get_sample_data
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
# -----------------------------------------------------------------------------

#SUBROUTINE-----------------------------------------------
def fermiDirac (mu, eig):
  kb = 8.617343 * 10**(-5)  # eV/K
  T  = 300.0   # K
  occ = 1.0 / (1.0 + (math.exp((eig - mu) / (kb * T))))
  return occ

def polygon_under_graph(x, y):
    """
    Construct the vertex list which defines the polygon filling the space under
    the (x, y) line graph. This assumes x is in ascending order.
    """
    return [(x[0], 0.), *zip(x, y), (x[-1], 0.)]

# New ----------------------------------------------------------

def new_wf_projection(figType,zbox,mu,muG,region,eigfile):
    global wfdata
    global verts, fcolor_list
    global Eig_list, Z_list, S_list, S_list_w, Occ_list, Occ_list_z, Eig_list_z
    global Left_Last_i
    global N_color
    global jet_colors

    for wi, a_wf in enumerate(wfdata):
        a_wf = str(a_wf)
        #print ('a_wf', a_wf)
        p_file = open(a_wf, 'r')
        rdata = p_file.readlines()
        p_file.close()
        # find the starting point of values
        startset = []
        spt = 664                                     # line where to start DATA in *.cube
        # x-axis number of grid
        xgrid = rdata[3]
        xgrid = xgrid.split()
        xgrid = int(xgrid[0])
        # y-axis number of grid
        ygrid = rdata[4]
        ygrid = ygrid.split()
        ygrid = int(ygrid[0])
        # z-axis number of grid
        zgrid = rdata[5]
        zgrid = zgrid.split()
        zgrid = int(zgrid[0])
        # eigen value
        eig = rdata[1]
        eig = eig.split('.')
        wfn = eig[2]
        wfn = wfn[2:]
        eig = eigfile[int(wfn)-1] 
        
        Eig_list.append(float(eig)-muG ) # eigen value list (zs value in plot)

        if region == 'left':            
            occ = fermiDirac(mu, float(eig))
            Left_Last_i = wi
        elif region == 'right':            
            occ = 1-fermiDirac(mu, float(eig))


        Occ_list.append(occ)
        #Eig_listreal.append(float(eig)-muG)
        
        print (" Eigenvalue, WFN = ", float(eig)+4.535625, " ", wfn)
        dummy = float(eig)+4.535625
        if region == 'left': 
            os.system('cp -rf %s %s+%s' %(a_wf, a_wf, dummy))
        elif region == 'right':
            os.system('cp -rf %s %s+%s' %(a_wf, a_wf, dummy))
# xy-plane summation routine
        value = rdata[spt:]   # wf data
        vvalue = []
        for b in value:
            for x in b.split():
                vvalue.append(x)
        zarray = []           # z-axis array
        svalue = []           # summation of WF along z-axis
        svalue_w = []
        eig_z  =[]
        occ_z  =[]

        disp = zbox/zgrid
        #tot = 0.0
        #for x in range(zgrid):
        #    for y in range(xgrid*ygrid):
        #        tot += abs(float(vvalue[x+zgrid*y]))**2

        Edge = 35

        for zi in range(zgrid):
            xr = zi*disp
            zarray.append(xr)

            eig_z.append(float(eig)-muG)
            occ_z.append(np.abs(occ))
            zv = 0.0             

            #if (zi >= 0 and zi<=Edge) or (zi <= zgrid-1 and zi >= zgrid-Edge) or (occ < 1e-11) or occ > (1-1e-11):
            if (zi >= 0 and zi<=Edge) or (zi <= zgrid-1 and zi >= zgrid-Edge):
            
                svalue.append(zv)
                continue

            for y in range(xgrid*ygrid):
                zv += abs(float(vvalue[zi+zgrid*y]))**2
            
            svalue.append(zv)
            svalue_w.append(zv*np.abs(occ) )
            
        
        Z_list.append(zarray)
        S_list.append(svalue)
        S_list_w.append(svalue_w)

        Eig_list_z.append(eig_z)

        Occ_list_z.append(occ_z)
            #verts.append(list(zip(zarray, svalue)))
        
# New ver 2 by JLee (190319)
        
    
    if region == 'left':         
        for a_occ in Occ_list[:Left_Last_i+1]:
            col_B = 0.0 
            col_R = 1.0 #float(a_occ) #1.0

            #if a_occ<=0.0 : col_G = 0.0
            #else          : col_G =   - np.log10(a_occ)/10

            #col_G = np.log10(a_occ+1)/np.log10(2)

            col_G =   - np.log10(a_occ)/10
            
            if col_G <= 0.0   : col_G = 0.0
            elif col_G >= 1.0 : col_G = 1.0

            a_col = (col_R, col_G, col_B, 1.0) #float(a_occ)+0.5)
            fcolor_list.append(a_col)
        
    else:        
        for a_occ in Occ_list[Left_Last_i+1:]:
            col_R = 0.0 #0.5 - float(a_occ)/2 
            col_B = 1.0 #float(a_occ) #1.0

            #if a_occ<=0.0 : col_G = 0.0
            #else          : col_G = -np.log10(a_occ)/10
            col_G =   - np.log10(a_occ)/10

            #col_G = np.log10(a_occ+1)/np.log10(2)
            
            if col_G <= 0.0: col_G = 0.0
            a_col = (col_R, col_G, col_B, 1.0) #float(a_occ)+0.5)
            fcolor_list.append(a_col)

        
# 사용자 정의 색상 함수
def color_function(xi, yi, zi, oi, mode):

    global N_color
    global jet_colors
    '''
    if   mode =='L': cut = 0.1
    elif mode =='R': cut = 0.1

    #if zi==0.0:
    #    ci = 0
    #else:    
    #    ci = np.log10(zi+1)/np.log10(2)
    #    #ci = -np.log10(zi+1)/10

    

    ci = (zi - cut)/(1 - cut)

    jet_i = int(ci * (N_color-1))
    jet_colors_1 = jet_colors[jet_i]

    c_r = jet_colors_1[0]
    c_g = jet_colors_1[1]
    c_b = jet_colors_1[2]

    if   mode =='L': 
        c_r = 1
        #c_b = 0
    elif mode =='R': 
        #c_r = 0
        c_b = 1
    '''
    alpha = 1

    c_r = 0
    c_g = 0
    c_b = 0

    F_si = 0.1


    if   mode =='L': 
        if zi >= F_si and zi< F_si + 0.1 :        
            c_r = 1
            c_g = zi
            c_b = 0
    elif mode =='R': 
        if zi >= F_si and zi< F_si + 0.1 :        
            c_r = 0
            c_g = zi
            c_b = 1
        

    return (c_r, c_g, c_b, alpha)

def x_scale(x, pos):    
    global x_max, x_min    
    global center_x_norm

    return f"{int((x - 0)*(x_max - x_min) + x_min )}"

def y_scale(y, pos):    

    global y_max, y_min
    

    return f"{round(y*(y_max - y_min) + y_min, 2)}"


def draw_QFL(dir_Loc):
    
    global wfdata

    global verts, fcolor_list
    global Eig_list, Z_list, S_list, S_list_w, Occ_list, Occ_list_z, Eig_list_z
    global Left_Last_i
    global N_color
    global jet_colors

    global x_max, x_min
    global y_max, y_min

    cwd = os.getcwd()
    print('draw_QFL-------')
    print (" Reading cube file .... \n")
    if not dir_Loc.endswith('/'):
        dir_Loc += '/'
    #### INPUT #####

    verts = []
    figType = 1   # type 1: 3D plot
                # type 2: 2D plot
    Eig_list_z = []

    Eig_list = []
    Z_list   = []
    S_list   = []
    S_list_w   = []
    Occ_list  = []
    Occ_list_z  = []
    fcolor_list = []

    Left_Last_i = 335

    zbox = 648.47900     # Ang

    #### INPUT(Chem. Pot) #####
    muL = -4.03390071818834         #float(sys.argv[3])
    muR = -5.03348199975444         #float(sys.argv[4])
    muG = -4.535625          #float(sys.argv[2])



    eigfilename = dir_Loc + "pn_junction.EIG"
    eigf = open(eigfilename, 'r')
    line = eigf.read()
    eigfile = line.split()

    N_wf = 700 # Wave function is 700
    N_z  = 500 #    Z position is 500

    #root_dir = os.getcwd()    
    left_right = ['left', 'right']

    for ER_i in left_right:
        #work_dir = root_dir + "/QFL_data/"+ER_i
        print(f"dir_Loc {dir_Loc}")
        work_dir = dir_Loc + ER_i
        os.chdir(work_dir)        

        wfdata = glob.glob("*.REAL.cube")   # read only real cube         
        wfdata = np.sort(wfdata)
        if ER_i == 'left':
            print ('JLee&HYeo: Left')
            mu = muL
            n(figType, zbox, mu,muG, ER_i, eigfile)
        else:
            print ('JLee&HYeo: Right')
            mu = muR
            n(figType, zbox, mu,muG, ER_i, eigfile)
        os.chdir(cwd)

    eigf.close()
    print (" plot generating .... \n")


    # X좌표는 z가 되어야 하고
    # Y좌표는 에너지가 되어야 한다. 
    
    Z_values_L = np.array(S_list[:Left_Last_i+1])  # 임의의 z 값 생성
    Z_values_R = np.array(S_list[Left_Last_i+1:])  # 임의의 z 값 생성
    
    Eig_list_L = Eig_list[:Left_Last_i+1]
    Eig_list_R = Eig_list[Left_Last_i+1:]

    #----------
    fig = plt.figure(figsize=(10, 10))
    ax  = fig.add_subplot(111, projection='3d', proj_type='ortho')

    # 축 레이블 설정
    ax.set_xlabel('$z(\AA)$')
    ax.set_ylabel('Energy(eV)')
    ax.set_zlabel('')

    #--------------------------------------------------------------

    Px = np.array(Z_list[0])


    Half_N_wf = int(N_wf/2)



    verts_L = [ polygon_under_graph(Px, Z_values_L[wi] ) for wi in range(Left_Last_i+1) ]
    verts_R = [ polygon_under_graph(Px, Z_values_R[wi] ) for wi in range(N_wf-Left_Last_i-1) ]

    facecolors_L = []
    for wi in range(Left_Last_i+1): facecolors_L.append(fcolor_list[wi])

    facecolors_R = []
    for wi in range(N_wf-Left_Last_i-1): facecolors_R.append(fcolor_list[wi+Left_Last_i+1])

    poly_1 = PolyCollection(verts_L, facecolors = facecolors_L, alpha=1.0) # edgecolors='white'
    poly_2 = PolyCollection(verts_R, facecolors = facecolors_R, alpha=1.0) # edgecolors='white'

    ax.add_collection3d(poly_1, zs=Eig_list_L, zdir='y')
    ax.add_collection3d(poly_2, zs=Eig_list_R, zdir='y')

    ax.set(xlim=(Px[0], Px[-1]), ylim=(Eig_list_L[0], Eig_list_R[-1]), zlim=(0, 0.06) )

    # Occ plot -----------------------------------------------------------
    x = np.linspace(Eig_list_L[0], Eig_list_R[-1], len(Occ_list[:Left_Last_i+1]))
    y = np.array(Occ_list[:Left_Last_i+1])*0.05
    #ax.plot(x, y, zs=0, zdir='x', label='curve in (x, y)')

    x = np.linspace(Eig_list_L[0], Eig_list_R[-1], len(Occ_list[Left_Last_i+1:]))
    y = np.array(Occ_list[Left_Last_i+1:])*0.05
    #ax.plot(x, y, zs=Px[-1], zdir='x', label='curve in (x, y)')
    # Occ plot end -----------------------------------------------------------

    ax.view_init(elev=91, azim=270)
    #ax.view_init(elev=5, azim=91)

    #work_dir = root_dir + "/QFL_data"
    #os.chdir(work_dir)
    #os.chdir(dir_Loc)

    fig_name= 'QFL_figure'

    plt.savefig(f'{fig_name}.png')

    plt.show()

    #plt.close()

    X_point, Y_point = np.meshgrid(Px, Eig_list)
    Z_value = np.array(S_list)  

    Z_value_max = np.max(Z_value)
    Z_value_min = np.min(Z_value)

    Z_value_norm = np.abs(Z_value/(Z_value_max-Z_value_min))

    Z_values_L_norm = Z_value_norm[:Left_Last_i+1]  # 임의의 z 값 생성
    Z_values_R_norm = Z_value_norm[Left_Last_i+1:]  # 임의의 z 값 생성

    Z_value_log = np.log10(Z_value_norm + 1)/np.log10(2)

    Z_value_log_L = Z_value_log[:Left_Last_i+1]
    Z_value_log_R = Z_value_log[Left_Last_i+1:]

    fig = plt.figure(figsize=(10, 10))

    ax  = fig.add_subplot(111, projection='3d', proj_type='ortho')

    # 축 레이블 설정

    ax.set_xlabel('$z(\AA)$')
    ax.set_ylabel('Energy(eV)')
    ax.set_zlabel('')

    surf = ax.plot_surface(X_point, Y_point, Z_value_norm, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.view_init(elev=90, azim=270)
    plt.show()
    #plt.close()


    # -------------------------------------------------------------------------------------------\

    N_color = 256

    jet = cm.get_cmap('jet', N_color)  
    jet_colors = jet(range(N_color))


    Px = np.array(Z_list[0])

    print('center -------------------------')

    center_x = (Px[Left_Last_i] + Px[Left_Last_i+1])/2
    print(center_x)


    x_max = np.array(Z_list).max()
    x_min = np.array(Z_list).min()

    xp_L = np.array(Z_list[:Left_Last_i+1]).flatten()  
    xp_L_norm = (xp_L - x_min)/(x_max - x_min) 

    xp_R = np.array(Z_list[Left_Last_i+1:]).flatten()  
    xp_R_norm = (xp_R - x_min)/(x_max - x_min) 

    y_max = np.array(Eig_list_z).max()
    y_min = np.array(Eig_list_z).min()

    yp_L = np.array(Eig_list_z[:Left_Last_i+1]).flatten()  
    yp_L_norm = (yp_L - y_min)/(y_max - y_min) 

    yp_R = np.array(Eig_list_z[Left_Last_i+1:]).flatten()
    yp_R_norm = (yp_R - y_min)/(y_max - y_min) 

    center_x_norm = (center_x - x_min)/(x_max - x_min) 

    print('center_x_norm ', center_x_norm)

    Z_value = np.array(S_list)  

    Z_value_max = np.max(Z_value)
    Z_value_min = np.min(Z_value)

    Z_value_norm = (Z_value-Z_value_min)/(Z_value_max-Z_value_min)

    Z_value_norm = np.log10(Z_value_norm + 1)/np.log10(2)

    Z_values_L_norm = Z_value_norm[:Left_Last_i+1].flatten()  # 임의의 z 값 생성
    Z_values_R_norm = Z_value_norm[Left_Last_i+1:].flatten()  # 임의의 z 값 생성

    print('------------------------------')

    Len_L = len(Z_values_L_norm)
    Len_R = len(Z_values_R_norm)

    L_max =  Z_values_L_norm.max()
    L_min =  Z_values_L_norm.min()

    R_max =  Z_values_R_norm.max()
    R_min =  Z_values_R_norm.min()


    occ_value = np.array(Occ_list_z)  

    occ_value_max = np.max(occ_value)
    occ_value_min = np.min(occ_value)

    occ_value_norm = (occ_value-occ_value_min)/(occ_value_max-occ_value_min)

    occ_values_L_norm = occ_value_norm[:Left_Last_i+1].flatten()  # 임의의 z 값 생성
    occ_values_R_norm = occ_value_norm[Left_Last_i+1:].flatten()  # 임의의 z 값 생성


    fig = plt.figure(figsize=(10, 10))

    ax  = fig.add_subplot(111)
    #fig.patch.set_facecolor('black')  # Figure의 배경색 설정
    ax.set_facecolor('black')  # 원하는 색으로 변경하세요.
    # 색상 계산
    colors_L = np.array([color_function(xi, yi, zi, oi, 'L') for xi, yi, zi, oi in zip(xp_L_norm, yp_L_norm, Z_values_L_norm, occ_values_L_norm)])
    colors_R = np.array([color_function(xi, yi, zi, oi, 'R') for xi, yi, zi, oi in zip(xp_R_norm, yp_R_norm, Z_values_R_norm, occ_values_R_norm)])

    # Scatter plot with x_scale markers
    ax.scatter(xp_L_norm, yp_L_norm, color=colors_L, marker='s', s=3, edgecolors='none')  # 's'는 사각형, s는 마커 크기
    ax.scatter(xp_R_norm, yp_R_norm, color=colors_R, marker='s', s=3, edgecolors='none')  # 's'는 사각형, s는 마커 크기

    ax.set_xlabel('$z(\AA)$')
    ax.set_ylabel('Energy(eV)')

    # FuncFormatter를 사용하여 함수를 눈금 포맷터로 설정합니다.
    x_formatter = FuncFormatter(x_scale)
    ax.xaxis.set_major_formatter(x_formatter)


    y_formatter = FuncFormatter(y_scale)
    ax.yaxis.set_major_formatter(y_formatter)

    ax.set(xlim=(0, 1))
    ax.set(ylim=(0, 1))

    plt.show()
    #plt.close()


    print('This is in oorinano')

