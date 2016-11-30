#!/usr/bin/env Python3.5

#!/usr/bin/env Python3.5

# Problem 2 (source term)

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize, xmax, ymax):
    afname = "sourc_a_" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" + str(ymax) + ".txt"
    fname = "sourc_c_" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" + str(ymax) + ".txt"

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

    meshsizes = [10,20,40,80,160,320,640]

    fs = 20.
    fw = "bold"
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    l2_norms = []
    
    xpos,se,asol,csol=read_files(meshsizes[0],xmax,ymax)
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(csol),np.max(csol),20)
    CS = ax1.contour(cols_plot, rows_plot, csol, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot, rows_plot, csol, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Source Terms for the Energy Equation, Meshsize = {}'.format(meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$FI_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    fig_name = 'mech510_p3_5.png'
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
    ax2.set_title('Contour Plot of the Analytic Source Terms for the Energy Equation, Meshsize = {}'.format(meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$FI_{analytic}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    fig_name = 'mech510_p3_6.png'
    fig2.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(se),np.max(se),20)
    CS = ax3.contour(cols_plot, rows_plot, se, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(cols_plot, rows_plot, se, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title('Contour Plot of the Source Term Error for the Energy Equation, Meshsize = {}'.format(meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$FI_{analytic} - FI_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    fig_name = 'mech510_p3_7.png'
    fig3.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    l2_norms.append(calc_l2_norm(se, meshsizes[0]))
    
    
    for meshsize in meshsizes[1:]:
        xpos,se,asol,csol=read_files(meshsize,xmax,ymax)
        l2_norms.append(calc_l2_norm(se, meshsize))
        
    
    fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
    
    l2_norms.sort()
    l2_norms = l2_norms[::-1]
    
    log_mesh = np.log10(meshsizes)
    log_norms = np.log10(l2_norms)
    
    reg = np.polyfit(log_mesh, log_norms, 1)
    fit = np.polyval(reg, log_mesh)
    
    ax4.plot(log_mesh,log_norms,'or')
    ax4.plot(log_mesh,fit,'-b')
    
    print(reg)
    
    ax4.set_xlabel(r'log$_{10}$(Meshsize)', fontsize=fs, fontweight=fw)
    ax4.set_ylabel(r'log$_{10}$(L$_2$-Norm)', fontsize=fs, fontweight=fw)
    ax4.set_title(r'Log-Log Plot of Source L$_2$-Norms for Different Mesh Sizes', fontsize=fs, fontweight=fw)
    ax4.tick_params(axis='both',which='major',labelsize=20)
    ax4.set_xlim([0.8,3.0])
    fig_name = 'mech510_p3_8.png'
    
    plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
       
    """
    fig0,ax0 = plt.subplots(1,1,figsize=(16,12))
    fig0a,ax0a = plt.subplots(1,1,figsize=(16,12))
    #fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    
    l2_norms = []
    for meshsize in meshsizes:
        xpos,sol_error,asolution,csolution=read_files(meshsize, xmax, tmax, cfl, ic)
        l2_norms.append(calc_l2_norm(sol_error,meshsize))
        
        if meshsize == 20:
            ax0.plot(xpos, csolution, label='Computed')
            ax0.plot(xpos, asolution, label='Analytic')
        
        if meshsize == 40:
            ax0a.plot(xpos, csolution, label='Computed')
            ax0a.plot(xpos, asolution, label='Analytic')
            
        if meshsize == 20 or meshsize == 40:
            label_str = "Meshsize = {}".format(meshsize)
            ax1.plot(xpos, sol_error, label=label_str)
    

    
    ax0.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax0.set_ylabel(r'$\overline{T_i}$', fontsize=fs, fontweight=fw)
    ax0.set_title("Computed and Analytic Solutions for Meshsize = 10, t = 1", fontsize=fs, fontweight=fw)    
    ax0.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p3_0.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    ax0a.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax0a.set_ylabel(r'$\overline{T_i}$', fontsize=fs, fontweight=fw)
    ax0a.set_title("Computed and Analytic Solutions for Meshsize = 40, t = 1", fontsize=fs, fontweight=fw)    
    ax0a.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p3_0a.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    """
        
    """
    ax1.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('Error', fontsize=fs, fontweight=fw)
    ax1.set_title("Solution Error (Analytic - Computed), t = 1", fontsize=fs, fontweight=fw)    
    ax1.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p2_1.png'
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    """

    """   
    l2_norms.sort()
    l2_norms = l2_norms[::-1]
    
    print(l2_norms)
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    
    log_mesh = np.log10(meshsizes)
    log_norms = np.log10(l2_norms)
    
    reg = np.polyfit(log_mesh, log_norms, 1)
    fit = np.polyval(reg, log_mesh)
    
    ax2.plot(log_mesh,log_norms,'or')
    ax2.plot(log_mesh,fit,'-b')
    
    print(reg)
    
    ax2.set_xlabel(r'log$_{10}$(Meshsize)', fontsize=fs, fontweight=fw)
    ax2.set_ylabel(r'log$_{10}$(L$_2$-Norm)', fontsize=fs, fontweight=fw)
    ax2.set_title(r'Log-Log Plot of L$_2$-Norms for Different Mesh Sizes', fontsize=fs, fontweight=fw)
    ax2.set_xlim([0.8,2.4])
    fig_name = 'mech510_p2_2.png'
    
    plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    print(10**((-3.0-reg[1])/reg[0]))
    print(10**((-4.0-reg[1])/reg[0]))
    """