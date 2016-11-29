#!/usr/bin/env Python3.5

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize, xmax, tmax, cfl, ic):
    afname = "analytic_" + str(meshsize) + "_xmax_" + str(xmax) + "_tmax_" + str(tmax) + "_cfl_" + str(cfl) + "_ic_" + str(ic) + ".txt"
    fname = "compute_" + str(meshsize) + "_xmax_" + str(xmax) + "_tmax_" + str(tmax) + "_cfl_" + str(cfl) + "_ic_" + str(ic) + ".txt"

    asolution = []
    with open(afname, 'r') as afile:
        for line in afile:
            asolution.append(float(line.split()[0]))
            
    asolution = np.array(asolution)
    
    csolution = []
    with open(fname, 'r') as cfile:
        for line in cfile:
            csolution.append(float(line.split()[0]))
            
    csolution = np.array(csolution)
    sol_error = asolution - csolution
    
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    return(xpos,sol_error,asolution,csolution)
    
    

def calc_l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/meshsize)
    
    return(l2_norm)

if __name__ == '__main__':
     
    xmax = 1
    tmax = 1
    cfl = int(0.4*100)
    ic = 1
    
    meshsizes = [10,20,40,80,160]

    fs = 20.
    fw = "bold"
   
    plt.close('all')
        
    plt.style.use('ggplot')
    
    fig0,ax0 = plt.subplots(1,1,figsize=(16,12))
    fig0a,ax0a = plt.subplots(1,1,figsize=(16,12))
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    
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
    ax0.set_title("Computed and Analytic Solutions for Meshsize = 20, t = 1", fontsize=fs, fontweight=fw)    
    ax0.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p2_0.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    ax0a.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax0a.set_ylabel(r'$\overline{T_i}$', fontsize=fs, fontweight=fw)
    ax0a.set_title("Computed and Analytic Solutions for Meshsize = 40, t = 1", fontsize=fs, fontweight=fw)    
    ax0a.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p2_0a.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    ax1.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('Error', fontsize=fs, fontweight=fw)
    ax1.set_title("Solution Error (Analytic - Computed), t = 1", fontsize=fs, fontweight=fw)    
    ax1.legend(loc='upper left',fontsize=fs)
    
    fig_name = 'mech510_p2_1.png'
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    

        
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
    
    
    plt.show()
    