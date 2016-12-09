# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 21:40:10 2016

@author: Timothy
"""

#!/usr/bin/env Python3.5

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize, xmax, tmax, cfl, ic, bc):
    afname = "analytic_" + str(meshsize) + "_xmax_" + str(xmax) + "_tmax_" + str(tmax) + "_cfl_" + str(cfl) + "_ic_" + str(ic) + "_bc_" + str(bc) + ".txt"
    fname = "compute_" + str(meshsize) + "_xmax_" + str(xmax) + "_tmax_" + str(tmax) + "_cfl_" + str(cfl) + "_ic_" + str(ic) + "_bc_" + str(bc) + ".txt"

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
    
    print("Opening {}".format(fname))
    
    return(xpos,sol_error,asolution,csolution)
    
    

def calc_l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/meshsize)
    
    return(l2_norm)
    
def calc_l1_norm(error, meshsize):
    error_sum = np.sum(np.abs(error))
    l1_norm = (error_sum/meshsize)
    
    return(l1_norm)
    
def calc_linf_norm(error,meshsize):
    abs_error = np.abs(error)
    
    linf_norm = np.max(abs_error)
    
    return(linf_norm)

if __name__ == '__main__':
     
    xmax = 1
    tmax = 1
    cfl = int(0.4*1000)
    ic = 1
    bcs = [1,2,3]
    
    meshsizes = np.array([10,20,40,80,160,320,640])
    logmesh = np.log10(meshsizes)

    fs = 20.
    fw = "bold"
   
    plt.close('all')
        
    plt.style.use('ggplot')
    
    fig0,ax0 = plt.subplots(1,1,figsize=(16,12))
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
 
    
    meshsize = meshsizes[3]
    
    bc_strs=["2nd-Upwind","2nd-Centred","1st-Upwind"]
    bcind = 0
    colors = ['c','g','r']
    for bc in bcs:
        l1_norms = []
        linf_norms = []
        l2_norms = []
        
        bc_str = bc_strs[bcind]
        color = colors[bcind]
        bcind = bcind + 1
        
        label1 = bc_str + r", L$_1$-Norm"
        labelinf = bc_str + r", L$_\infty$-Norm"
        for meshsize in meshsizes:
            xpos,sol_error,asolution,csolution=read_files(meshsize, xmax, tmax, cfl, ic ,bc)
            l1_norms.append(calc_l1_norm(sol_error,meshsize))
            linf_norms.append(calc_linf_norm(sol_error,meshsize))
            l2_norms.append(calc_l2_norm(sol_error,meshsize))
            
            if meshsize == 80:
                ax0.plot(xpos, csolution, color=color, label=bc_str)
                ax1.plot(xpos, sol_error, color=color, label=bc_str)
                
        logl1 = np.log10(np.array(l1_norms))
        loglinf = np.log10(np.array(linf_norms))
        logl2 = np.log10(np.array(l2_norms))
        
        ax2.plot(logmesh,logl1,'o',markersize=10,color=color)
        ax2.plot(logmesh,loglinf,'x',markersize=10,color=color)
        
        reg = np.polyfit(logmesh, logl1, 1)
        fit = np.polyval(reg, logmesh)        
        ax2.plot(logmesh,fit,'-',color=color, label=label1)
        
        print(bc)
        print(reg)
            
        reg = np.polyfit(logmesh, loglinf, 1)
        fit = np.polyval(reg, logmesh)        
        ax2.plot(logmesh,fit,'--',color=color, label=labelinf)
        
        reg = np.polyfit(logmesh, logl2, 1)
        print(bc)
        print(reg)

 
    ax0.set_xlabel(r'$x$', fontsize=fs+20, fontweight=fw)
    ax0.set_ylabel(r'$\overline{T_i}$', fontsize=fs+20, fontweight=fw)
    ax0.set_title("Computed Solutions with Different Flux Calculations at i = 3/2, Meshsize = 80", fontsize=fs, fontweight=fw)    
    ax0.legend(loc='upper left',fontsize=fs)
    ax0.tick_params(axis='both',which='major',labelsize=20)
    
    fig_name = 'mech510_p2_5.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    ax1.set_xlabel(r'$x$', fontsize=fs+20, fontweight=fw)
    ax1.set_ylabel('Error', fontsize=fs, fontweight=fw)
    ax1.set_title("Solution Error with Different Flux Calculations at i = 3/2, Meshsize = 80", fontsize=fs, fontweight=fw)    
    ax1.legend(loc='upper left',fontsize=fs)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    
    fig_name = 'mech510_p2_6.png'
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
        
    ax2.set_xlabel(r'log$_{10}$(Meshsize)', fontsize=fs, fontweight=fw)
    ax2.set_ylabel(r'log$_{10}$(Error Norm)', fontsize=fs, fontweight=fw)
    ax2.set_title(r'Log-Log Plot of Error Norms for Different Mesh Sizes and Boundary Conditions', fontsize=fs, fontweight=fw)
    ax2.set_xlim([0.8,3.0])
    handles, labels = ax2.get_legend_handles_labels()
    new_labels = [labels[5],labels[1],labels[3],labels[4],labels[0],labels[2]]
    new_handles = [handles[5],handles[1],handles[3],handles[4],handles[0],handles[2]]
    ax2.legend(new_handles,new_labels,loc='upper right',fontsize=fs)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    
    
    fig_name = 'mech510_p2_7.png'
    fig2.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)