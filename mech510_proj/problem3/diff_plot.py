#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import sys


def savefigure(fig_name,fignum):
    print("Saving to {}".format(fig_name))
    fignum.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    return(None)
    
    
    
if __name__ == '__main__':
    
        
    #plt.close('all')
    plt.style.use('ggplot')
    
    timings = []
    

    fs = 20.
    fw = "bold"

    upfilename = "u100.txt"
    umfilename = "u_100.txt"

    meshsize = 20
    cols = np.linspace(0,1,meshsize+1)
    rows = np.linspace(0,1,meshsize+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = -(rows + np.diff(rows)[0]/2)[:-1]

    uparray = []
    umarray = []

    with open(upfilename, 'r') as ufile:
        for line in ufile:
            num_list = [float(num) for num in line.split()]
            uparray.append(num_list)    

    
    with open(umfilename, 'r') as ufile:
        for line in ufile:
            num_list = [float(num) for num in line.split()]
            umarray.append(num_list)    

    uparray = np.flipud(np.array(uparray))

    
    umarray = np.fliplr(np.flipud(np.array(umarray)))
    
    contours = uparray + umarray
    
    
    # u-velocity plot
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(contours),np.max(contours),80)
    contour_filled = ax2.contourf(cols_plot, rows_plot, contours, levels)
    CS = ax2.contour(cols_plot, rows_plot, contours, levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title(r'u(x,y) with Tw = 100 Plus u(1-x,y) with Tw = -100',fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'u(x,y,T$_w$ = 100) + u(1-x,y,T$_w$ = -100)',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    figname="diff_u.png"
    savefigure(figname,fig2)
    
    plt.show()
    

   
