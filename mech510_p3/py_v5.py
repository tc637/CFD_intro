#!/usr/bin/env Python3.5

#!/usr/bin/env Python3.5

# Problem 5 

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize, xmax, ymax):
    fname = ("Tmesh_cn" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" 
             + str(ymax) + "_t_" + "0" + "_scheme_" + "2" + ".txt")
    afname = ("Tmesh_cn" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_"  
             + str(ymax) + "_t_" + "1" + "_scheme_" + "2" + ".txt")

    asolution = []
    with open(afname, 'r') as afile:
        for line in afile:
            num_list = [float(num) for num in line.split()]
            asolution.append(num_list)
            
    asolution = np.array(asolution)
    
    csolution = []
    with open(fname, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            csolution.append(num_list)
            
    csolution = np.array(csolution)
    sol_error = asolution - csolution
    
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    return(xpos,sol_error,asolution,csolution)
    
    

def calc_l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/(meshsize*meshsize))

    return(l2_norm)
    
    
def read_analytic(thefile, asolution):
    
    with open(thefile, 'r') as afile:
        for line in afile:
            line = line.split()
            if len(line) == 0:
                continue
            else:
                numline = []
                for each_part in line:
                    numline.append(float(each_part))
                    
                asolution[int(numline[1]),int(numline[0])] = numline[2]
    
    return(None)

if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
     
    xmax = 1
    ymax = 1

    meshsizes = [10,20,40,80,160,320,640]
    ind = 0

    fs = 20.
    fw = "bold"
    
    cols = np.linspace(0,1,meshsizes[ind]+1)
    rows = np.linspace(0,1,meshsizes[ind]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    l2_norms = []
    
    ana = np.zeros([12,27])
    before = np.zeros([12,27])
    
    read_analytic("analytic.dat",ana)
    read_analytic("before.dat",before)
    
    
    cols = np.linspace(0,1,27)
    rows = np.linspace(0,1,12)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]*5
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    before = before[1:12,1:27]
    ana = ana[1:12,1:27]
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(before),np.max(before),20)
    CS = ax1.contour(cols_plot,rows_plot,before, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot,rows_plot,before, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    """
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(ana),np.max(ana),20)
    CS = ax2.contour(cols_plot,rows_plot,ana, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(cols_plot,rows_plot,ana, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    """
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(ana),np.max(ana),20)
    CS = ax2.contour(ana, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(ana, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    
   