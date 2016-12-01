#!/usr/bin/env Python3.5

#!/usr/bin/env Python3.5

# Problem 5 

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize,jmax,xmax,ymax,tmax,scheme,p):
    fname = ("Tmesh_cn" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" 
             + str(ymax) + "_t_" + str(tmax) + "_scheme_" + str(scheme) + "_p_" + str(p) + ".txt")

    csolution = []
    with open(fname, 'r') as cfile:
        for line in cfile:
            num_list = [float(num) for num in line.split()]
            csolution.append(num_list)
            
    csolution = np.array(csolution)
    
    dx = float(xmax)/float(meshsize)
    dy = float(ymax)/float(jmax)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    ypos = np.array([dy*(j+0.5) for j in range(0,jmax)])
    return(xpos,ypos,csolution)
    
    

def calc_l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/(meshsize*meshsize))

    return(l2_norm)
    
    
def read_analytic(thefile, asolution):
    
    with open(thefile, 'r') as afile:
        for line in afile:
            line = line.split()
            if len(line) == 0:
                continue
            else:
                numline = []
                for each_part in line:
                    numline.append(float(each_part))
                    
                asolution[int(numline[1]),int(numline[0])] = numline[2]
    
    return(None)

if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
     
    xmax = 5
    ymax = 1
    tmax = 0
    scheme = 2
    p = 2
    jmax = 10
    
    meshsizes = [25,200]
    ind = 0

    fs = 20.
    fw = "bold"
    
    cols = np.linspace(0,1,meshsizes[ind]+1)
    rows = np.linspace(0,1,meshsizes[ind]+1)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    l2_norms = []
    
    ana = np.zeros([12,27])
    before = np.zeros([12,27])
    
    read_analytic("analytic.dat",ana)
    read_analytic("before.dat",before)
    
    cols = np.linspace(0,1,26)
    rows = np.linspace(0,1,11)
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]*5
    rows_plot = (rows + np.diff(rows)[0]/2)[:-1]
    
    before = before[1:11,1:26]
    ana = ana[1:11,1:26]
    
    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(before),np.max(before),20)
    CS = ax1.contour(cols_plot,rows_plot,before, levels, colors='k')
    ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax1.contourf(cols_plot,rows_plot,before, levels)
    cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 

    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(ana),np.max(ana),20)
    CS = ax2.contour(cols_plot,rows_plot,ana, levels, colors='k')
    ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax2.contourf(cols_plot,rows_plot,ana, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    
    # ========================
    
    
    xpos,ypos,csol = read_files(meshsizes[ind],jmax,xmax,ymax,tmax,scheme,p)
    
    
    
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(csol),np.max(csol),20)
    CS = ax3.contour(xpos,ypos,csol, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(xpos,ypos,csol, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(meshsizes[ind]),
                  fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    
    ana_wall = ana[:,0]
    fig4,ax4=plt.subplots(1,1)
    ax4.plot(ana_wall)
    
    computed_wall = [0.20482879493674047,0.50441279748591161,0.69059679858669065,0.81781279893430447,0.92234879899224020,1.0223487989922411,1.1178127989343045,1.1905967985866905,1.2044127974859113,1.1048287949367397]
   
    ax4.plot(computed_wall)
    computed_wall = [0.21249274965409554,
  0.50905274923567090,
  0.69296874905703576,
  0.81867274900231868,
  0.92245274899427210,
   1.0224527489942721,
   1.1186727490023187,
   1.1929687490570358,
   1.2090527492356709,
   1.1124927496540955]
   

    ax4.plot(computed_wall)
    

   
   
   
   
   
   
   
  
  

    