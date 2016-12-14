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
    
    

if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
    xmax = 1
    ymax = 1

    meshsizes = [10]
    varind = [1,2,3]
    ind = 2
    
    fs = 20.
    fw = "bold"
    varind = 0
    
    varstrs = ['Pressure','u-Velocity','v-Velocity']

    varstr = varstrs[varind]
    
    cols = np.linspace(0,1,meshsizes[0]+1)
    rows = np.linspace(0,1,meshsizes[0]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    init_array = np.zeros([10,10,3])
    final_array = np.zeros([10,10,3]) # (j,i,variables)
    
    with open('test_data_init.txt','r') as init:
        for line in init:
            line = line.split()
            
            if len(line) == 0:
                continue
            else:
                iind = int(line[1])
                jind = int(line[3])
                
                if (jind == 0) or (jind == 11): # Exclude ghost cells
                    continue
                
                dP = float(line[5])
                dU = float(line[7])
                dV = float(line[9])
        
                init_array[jind-1,iind-1,0] = dP
                init_array[jind-1,iind-1,1] = dU
                init_array[jind-1,iind-1,2] = dV 
    
    with open('test_data_final.txt','r') as final:
        for line in final:
            line = line.split()
            
            if len(line) == 0:
                continue
            else:
                iind = int(line[1])
                jind = int(line[3])
                dP = float(line[5])
                dU = float(line[7])
                dV = float(line[9])
        
                final_array[jind-1,iind-1,0] = dP
                final_array[jind-1,iind-1,1] = dU
                final_array[jind-1,iind-1,2] = dV
             
    # Update final array
    final_array = final_array + init_array
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(init_array[:,:,varind]),np.max(init_array[:,:,varind]),20)
    CS = ax1.contour(cols_plot, rows_plot, init_array[:,:,varind], levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot, rows_plot, init_array[:,:,varind], levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Test Data for {}, t = 0, Meshsize = {}'.format(varstr,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{test}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(final_array[:,:,varind]),np.max(final_array[:,:,varind]),20)
    CS = ax2.contour(cols_plot, rows_plot, final_array[:,:,varind], levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(cols_plot, rows_plot, final_array[:,:,varind], levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Test Data for {}, t = 0.05, Meshsize = {}'.format(varstr,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{test}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    print(final_array[:,:,0])
    #print(final_array[4,4,1])
    #print(final_array[4,4,2])
    
    
    pfile0 = 'pmesh_cn10_xmax_1_ymax_1_t_0_p_2.txt'
    pfile1 = 'pmesh_cn10_xmax_1_ymax_1_t_1_p_2.txt'
    
    solution = []
    with open(pfile1, 'r') as pfile:
        for line in pfile:
            num_list = [float(num) for num in line.split()]
            solution.append(num_list)
            
    psolution = np.array(solution)
    
    
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(solution),np.max(solution),20)
    CS = ax3.contour(cols_plot, rows_plot, solution, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(cols_plot, rows_plot, solution, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title('Contour Plot of the Computed Solutions for {}, t = 0.05, Meshsize = {}'.format(varstr,meshsizes[0]),
                  fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$P_{test}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    

    plt.show()
    
    
    