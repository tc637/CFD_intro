#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize, xmax, ymax, varind):
    afname = "FI_anal_" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" + str(ymax) + "_var_" + str(varind) + ".txt"
    fname = "FI_comp_" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" + str(ymax) + "_var_" + str(varind) + ".txt"

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

if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
     
    xmax = 1
    ymax = 1

    meshsizes = [20, 40]
    varind = [1,2,3]
    ind = 2
    
    fs = 20.
    fw = "bold"
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    l2_norms = []
    
    xpos,se,asol,csol=read_files(meshsizes[0],xmax,ymax,varind[ind])
    
    varstrs = ['Pressure','u-Velocity','v-Velocity']
    varstr = varstrs[ind]
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(csol),np.max(csol),20)
    CS = ax1.contour(cols_plot, rows_plot, csol, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot, rows_plot, csol, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Computed {} Flux Integrals for the Navier-Stokes Equations, Meshsize = {}'.format(varstr,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$FI_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'mech510_proj_1c_' + str(varind[ind]) + ".png"
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(asol),np.max(asol),20)
    CS = ax2.contour(cols_plot, rows_plot, asol, levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(cols_plot, rows_plot, asol, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Analytic {} Flux Integrals for the Navier-Stokes Equations, Meshsize = {}'.format(varstr,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$FI_{analytic}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    fig_name = 'mech510_proj_1a_' + str(varind[ind]) + ".png"
    fig2.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    l2_norms.append(calc_l2_norm(se,meshsizes[0]))
    
    xpos,se,asol,csol=read_files(meshsizes[1],xmax,ymax,varind[ind])
    l2_norms.append(calc_l2_norm(se,meshsizes[1]))
    
    textfile = 'l2_norms_{}.txt'.format(varind[ind])
    with open(textfile,'w') as l2_file:
        for l2 in l2_norms:
            l2_file.write(str(l2)+'\n')
            print(l2)
        l2_file.write(str(l2_norms[1]/l2_norms[0])+'\n')    
        print(l2_norms[1]/l2_norms[0])
        
    #plt.show()
    
