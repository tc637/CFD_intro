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
    
    timings = []
    overs = np.arange(100,210,5)
    meshsizes=[10]
    fs = 20.
    fw = "bold"
    
    with open('overs.txt','r') as overfile:
        for line in overfile:
            line = line.split()
            
            if len(line) == 0:
                continue
            
            if line[0] == 'real':
                timing = line[1]
                timing = timing.split('m')
                timing = timing[1].split('s')[0]
                
                timings.append(np.float(timing))
    
    timings = np.array(timings)
    
    iterations = []
    pfiles = []
    ufiles = []
    vfiles = []
    
    
    with open('overs.stderr','r') as stdfile:
        for line in stdfile:
            line = line.split()

            if len(line) == 0:
                continue
            
            if line[0] == 'Writing':
                prefix = line[3].split('_')[0]
                
                if prefix == 'pmesh':
                    pfiles.append(line[3])
                elif prefix == 'umesh':
                    ufiles.append(line[3])
                elif prefix == 'vmesh':
                    vfiles.append(line[3])
                    
            if line[0] == 'Iteration':
                iterations.append(np.float(line[2]))
                
                
    iterations = np.array(iterations)
    pfiles = np.array(pfiles)
    ufiles = np.array(ufiles)
    vfiles = np.array(vfiles)
    
    itinds = np.where(iterations == np.min(iterations))
    timeinds = np.where(timings == np.min(timings))
    
    print(itinds)
    print(timeinds)
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    
    ax2 = ax1.twinx()

    ax1.plot(overs/100, iterations, '-b', label='Iteration Count')
    ax2.plot(overs/100, timings, '-r', label='Time to Steady State')
    
    ax1.set_xlabel(r'$\omega$', fontsize=fs+10, fontweight=fw)
    ax1.set_ylabel('Iteration Count', fontsize=fs, fontweight=fw, color='b')
    ax2.set_ylabel('Time', fontsize=fs, fontweight=fw, color='r')
    
    ax1.set_title(r'Iteration Count and Run Time with Different $\omega$',fontsize=fs,fontweight=fw)
    
    ax1.tick_params(axis='both',which='major',labelsize=20)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    
    fig_name = 'overs.png'
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    
    
    iterations = iterations[0:40]
    
    itinds = np.where(iterations == np.min(iterations))
    timeinds = np.where(timings == np.min(timings))
    
    print(itinds)
    print(timeinds)
    
    file1 = pfiles[0]
    print(file1)
    
        
    xmax = 1
    ymax = 1
    
    meshsize = meshsizes[0]
    
    solution = []
    with open(file1, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            solution.append(num_list)
            
    solution = np.array(solution)
    
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]

    
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),40)
    CS = ax3.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title(r'Contour Plot of the Steady-State Solution, $\omega$ = 1',fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'overs_ref.png'
    fig3.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
 
    ind = 0
    file1 = pfiles[itinds[0][ind]]
    print(file1)
    old_sol = solution
    solution = []
    with open(file1, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            solution.append(num_list)
            
    solution = np.array(solution)
    
    print(old_sol - solution)
    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]

    
    fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),40)
    CS = ax4.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax4.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax4.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax4.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax4.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax4.set_title(r'Contour Plot of the Steady-State Solution, $\omega$ = {}'.format(overs[itinds[0][ind]]/100),fontsize=fs,fontweight=fw)
    ax4.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{computed}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig_name = 'overs_contour.png'
    fig4.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    
    plt.show()
    
    
    
    
