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
    
    
    argus = sys.argv
    
    if len(sys.argv) < 10:
        print("meshsize,ymax,tlength,beta,apres,omega,Tw,dt")
        exit()
    
    meshsize = int(sys.argv[1])
    ymax = int(sys.argv[2])
    tmax = int(sys.argv[3])
    beta = int(sys.argv[4])
    apres = int(sys.argv[5])
    omega = int(sys.argv[6])
    Tw = int(sys.argv[7])
    dt = float(sys.argv[8])
    tol = int(sys.argv[9])
        
    #plt.close('all')
    plt.style.use('ggplot')
    
    timings = []
    

    fs = 20.
    fw = "bold"

    ufilename = "umesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"
    vfilename = "vmesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"
    pfilename = "pmesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"
    
    cols = np.linspace(0,1,meshsize+1)
    rows = np.linspace(0,1,meshsize+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = -(rows + np.diff(rows)[0]/2)[:-1]

    uarray = []
    varray = []
    parray = []
    
    with open(ufilename, 'r') as ufile:
        for line in ufile:
            num_list = [float(num) for num in line.split()]
            uarray.append(num_list)    
    with open(vfilename, 'r') as vfile:
        for line in vfile:
            num_list = [float(num) for num in line.split()]
            varray.append(num_list)
    with open(pfilename, 'r') as pfile:
        for line in pfile:
            num_list = [float(num) for num in line.split()]
            parray.append(num_list)

    uarray = np.flipud(np.array(uarray))
    varray = np.flipud(np.array(varray))
    speed = np.sqrt(uarray**2 + varray**2)
    
    parray = np.flipud(np.array(parray))
    
    
    # Pressure and velocity plot
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(parray),np.max(parray),80)
    CS = ax1.contour(cols_plot, rows_plot, parray, levels, colors='w',linestyles='dashed')
    #ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot, rows_plot, parray, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title(r'Steady-State Pressure Contours and Velocity Streamlines, Tw = {}, dt = {}'.format(Tw,dt),fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    X, Y = np.meshgrid(cols_plot,rows_plot)
    lw = 10*speed/speed.max() 
    Q = ax1.streamplot(X,Y,uarray,varray,color='k',linewidth=lw,zorder=5)
    
    figname = "pcontour_velocities_{}_{}_{}_{}_{}_{}_{}_{}.png".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize)
    savefigure(figname,fig1)
    
    # ===========================
    
    # u-velocity plot
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(uarray),np.max(uarray),80)
    contour_filled = ax2.contourf(cols_plot, rows_plot, uarray, levels)
    CS = ax2.contour(cols_plot, rows_plot, uarray, levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title(r'Steady-State u-Velocity Contours, Tw = {}, dt = {}'.format(Tw, dt),fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$u_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    figname = "ucontour_{}_{}_{}_{}_{}_{}_{}_{}.png".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize)
    savefigure(figname,fig2)
    
    # ===========================
    
    # v-velocity plot
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(varray),np.max(varray),80)
    contour_filled = ax3.contourf(cols_plot, rows_plot, varray, levels)
    CS = ax3.contour(cols_plot, rows_plot, varray, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title(r'Steady-State v-Velocity Contours, Tw = {}, dt = {}'.format(Tw, dt),fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$v_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    figname = "vcontour_{}_{}_{}_{}_{}_{}_{}_{}.png".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize)
    savefigure(figname,fig3)
    
    # ===========================
    
    uvels = uarray[:,int(np.shape(uarray)[1]/2)]
    
    fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
    ax4.plot(rows_plot, uvels)
    ax4.set_xlabel('y', fontsize=fs, fontweight=fw)
    ax4.set_ylabel('u', fontsize=fs, fontweight=fw)
    ax4.set_title(r'u-Velocity Along the Vertical Symmetry Line, Tw = {}'.format(Tw),fontsize=fs,fontweight=fw)
    ax4.tick_params(axis='both',which='major',labelsize=20)
    
    
    figname = "uslice_{}_{}_{}_{}_{}_{}_{}_{}.png".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize)
    savefigure(figname,fig4)
    
    # ===========================
    
    ufilename = "converg_"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+"_"+"1"+".txt"
    vfilename = "converg_"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+"_"+"2"+".txt"
    pfilename = "converg_"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+"_"+"3"+".txt"
    
    uconv= []
    vconv = []
    pconv = []
    
    with open(ufilename, 'r') as ufile:
        for line in ufile:
            uconv.append(np.float(line))    
    with open(vfilename, 'r') as vfile:
        for line in vfile:
            vconv.append(np.float(line))   
    with open(pfilename, 'r') as pfile:
        for line in pfile:
            pconv.append(np.float(line))   
      
    iterations = np.arange(1,len(uconv)+1)
    
    uconv = np.log10(np.array(uconv))
    vconv = np.log10(np.array(vconv))
    pconv = np.log10(np.array(pconv))
    
    fig5,ax5 = plt.subplots(1,1,figsize=(16,12))
    ax5.plot(iterations,pconv,'-b',label='P')
    ax5.plot(iterations,uconv,'-g',label='u')
    ax5.plot(iterations,vconv,'-r',label='v')
    
    ax5.set_xlabel('Iteration Count', fontsize=fs, fontweight=fw)
    ax5.set_ylabel(r'log$_{10}$(L$_{2}$-Norm)', fontsize=fs, fontweight=fw)
    ax5.set_title(r'Change in Solution vs Iteration Count, Tw = {}, dt = {}'.format(Tw, dt),fontsize=fs,fontweight=fw)
    ax5.tick_params(axis='both',which='major',labelsize=20)
    ax5.legend(fontsize=fs)
    
    figname = "converge_{}_{}_{}_{}_{}_{}_{}_{}.png".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize)
    savefigure(figname,fig5)
            
    plt.show()