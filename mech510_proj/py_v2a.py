#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt

def read_files(varmesh, meshsize, xmax, ymax, tmax):
    fname = varmesh + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" + str(ymax) + "_t_" + str(tmax) + "_p_" + str(2) + ".txt"
    print(fname)
    
    csolution = []
    with open(fname, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            csolution.append(num_list)
            
    csolution = np.array(csolution)
    
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    return(xpos,csolution)
    
    

if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
    xmax = 1
    ymax = 1

    meshsizes = [10]
    varind = [1,2,3]
    ind = 2
    
    fs = 20.
    fw = "bold"
    varind = 0
    
    varstrs = ['Pressure','u-Velocity','v-Velocity']
    varmeshes = ['pmesh_cn','umesh_cn','vmesh_cn']
    
    varstr = varstrs[varind]
    varmesh = varmeshes[varind]
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    timing = 31  ######
     
    xpos,solution=read_files(varmeshes[varind],meshsizes[0],1,1,timing)
    
    solution = np.array(solution)

    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),40)
    CS = ax1.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Computed Solutions for {}, t = {}, Meshsize = {}'.format(varstrs[varind],timing,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'mech510_proj_2converge_' + str(varind+1) + str(meshsizes[0]) +".png"
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    varind = 1
    
    xpos,solution=read_files(varmeshes[varind],meshsizes[0],1,1,timing)
    
    solution = np.array(solution)

    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),40)
    CS = ax2.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Computed Solutions for {}, t = {}, Meshsize = {}'.format(varstrs[varind],timing,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$u_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'mech510_proj_2converge_' + str(varind+1) + str(meshsizes[0]) +".png"
    fig2.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    varind = 2
    
    xpos,solution=read_files(varmeshes[varind],meshsizes[0],1,1,timing)
    
    solution = np.array(solution)

    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),40)
    CS = ax3.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title('Contour Plot of the Computed Solutions for {}, t = {}, Meshsize = {}'.format(varstrs[varind],timing,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$v_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'mech510_proj_2converge_' + str(varind+1) + str(meshsizes[0]) +".png"
    fig3.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    
    
    """
    fig_name = 'mech510_proj_2c_' + str(varind+1) + ".png"
    fig3.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    print(fig_name)
    """
    plt.show()
    
    
    