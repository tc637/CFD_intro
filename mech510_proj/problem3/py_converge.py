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
    
    #plt.close('all')
    plt.style.use('ggplot')
    
    timings = []
    overs = np.arange(100,200,5)
    meshsizes=[10]
    fs = 20.
    fw = "bold"
    
    meshsize = meshsizes[0]

    file1 = "converg_10_xmax_1_ymax_1_var_1.txt"
    file2 = "converg_10_xmax_1_ymax_1_var_2.txt"
    file3 = "converg_10_xmax_1_ymax_1_var_3.txt"
    
    dt = 0.5
    xmax = 1
    ymax = 1
    
    converge = []
    with open(file1, 'r') as cfile:
        for line in cfile:
            converge.append(np.float(line))

    fig1, ax1 = plt.subplots(1,1)
    
    ax1.plot(np.log10(converge), '-b', label='P')
        
    converge = []
    with open(file2, 'r') as cfile:
        for line in cfile:
            converge.append(np.float(line))
            
    ax1.plot(np.log10(converge), '-g', label='u') 
    
    
    converge = []   
    with open(file3, 'r') as cfile:
        for line in cfile:
            converge.append(np.float(line))
            
    ax1.plot(np.log10(converge), '-r', label='v')
    
    ax1.set_title(r'Change in Solution vs Iteration Count, dt = {}'.format(dt),fontsize=fs,fontweight=fw)
            
    ax1.set_xlabel('Iteration Count',fontsize=fs,fontweight=fw)
    ax1.set_ylabel(r'log$_{10}$(L$_{2}$-Norm)',fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
                            
    fig_name = 'l2_{}.png'.format(dt)
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    plt.legend(fontsize=fs)

    plt.show()
    
   
    
    
