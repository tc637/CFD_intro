#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage

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
    
    #plt.close('all')
    plt.style.use('ggplot')
    
    timings = []
    overs = np.arange(100,200,5)
    meshsizes=[41]
    fs = 20.
    fw = "bold"
    
    meshsize = meshsizes[0]

    file1 = "FILE1"
    dt = 0.5
    xmax = 1
    ymax = 1
    solution = []
    with open(file1, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            solution.append(num_list)
            
    solution = np.flipud(np.array(solution))
    
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    filt = 10
    
    print(np.shape(solution))
    solution = scipy.ndimage.zoom(solution, filt)
    print(np.shape(solution))
    
    cols = np.linspace(0,1,meshsizes[0]*filt+1)
    rows = np.linspace(0,1,meshsizes[0]*filt+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = -(rows + np.diff(rows)[0]/2)[:-1]
    
    fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),160)
    #CS = ax4.contour(cols_plot, rows_plot, solution, levels, colors='k')
    #ax4.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax4.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax4.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax4.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax4.set_title(r'Contour Plot of the Steady-State Solution for {}, dt = {}'.format("VAR",dt),fontsize=fs,fontweight=fw)
    ax4.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$VAR_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    plt.draw()
    """
    fig_name = 'overs_contour_{}.png'.format(float(OMEGA))
    fig4.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    """
    
    plt.show()
    
   
    
    
