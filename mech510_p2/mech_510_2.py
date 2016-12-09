#!/usr/bin/env Python3.5

import numpy as np
import matplotlib.pyplot as plt

def l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    
    l2_norm = np.sqrt(error_sum/meshsize)
    
    return(l2_norm)

if __name__ == '__main__':
    
    cfl = int(0.4*100)
    xmax = 1
    meshsize = 20

    fs = 20.
    fw = "bold"
    dx = xmax/meshsize
    
    plt.close('all')
        
    plt.style.use('ggplot')
    
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    afname = "analytic_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    fname = "compute_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    
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
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    
    label_str = "Meshsize = {}".format(meshsize)
    ax1.plot(xpos, sol_error, '-b', label=label_str)
    ax1.set_xlabel('t', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('Error', fontsize=fs, fontweight=fw)
    ax1.set_title("Time Series of the Solution Error (Analytical - Computed)", fontsize=fs, fontweight=fw)
    
    l2_norms = []
    l2_norms.append(l2_norm(sol_error,meshsize))
    
    cfl = int(0.4*100)
    xmax = 1
    meshsize = 40
    dx = xmax/meshsize

    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    afname = "analytic_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    fname = "compute_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    
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
    
    label_str = "Meshsize = {}".format(meshsize)
    ax1.plot(xpos, sol_error, '-r', label=label_str)
    
    ax1.legend(loc='upper right',fontsize=fs)
    
    fig_name = 'mech510_p2_1.png'
    plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
        
    l2_norms.append(l2_norm(sol_error,meshsize))
    
    xmax = 1
    meshsize = 10
    dx = xmax/meshsize
    
    afname = "analytic_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    fname = "compute_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    
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
    l2_norms.append(l2_norm(sol_error,meshsize))
    
    xmax = 1
    meshsize = 80
    dx = xmax/meshsize
    
    afname = "analytic_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    fname = "compute_" + str(meshsize) + "_max_" + str(xmax) + "_cfl_" + str(cfl) +".txt"
    
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
    l2_norms.append(l2_norm(sol_error,meshsize))
    
    meshsizes=[10,20,40,80]
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
    ax2.set_title(r'Log-Log Plot of L$_2$-Norms for Different Mesh Sizes for the Wave Equation', fontsize=fs, fontweight=fw)
    
    fig_name = 'mech510_p2_2.png'
    
    plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    plt.show()
    