#!/usr/bin/env Python3.5

#!/usr/bin/env Python3.5

# Problem 5 

import numpy as np
import matplotlib.pyplot as plt

def read_files(meshsize,jmax,xmax,ymax,tmax,scheme,p):
    fname = ("Tmesh_cn" + str(meshsize) + "_xmax_" + str(xmax) + "_ymax_" 
             + str(ymax) + "_t_" + str(tmax) + "_scheme_" + str(scheme) + "_p_" + str(p) + ".txt")

    print(fname)
    #fname="Tmesh_cn25_xmax_5_ymax_1_t_1_scheme_1_p_2.txt"
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
    
    

def calc_l2_norm(error, imax, jmax):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/(imax*jmax))

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

def analytic_sol(ypos):
    
    Pr = 0.7
    Ec = 0.1
    ubar = 3
    
    Tw = ypos + 3./4.*Pr*Ec*ubar**2*(1.-(1.-2.*ypos)**4)
    
    return(Tw)



if __name__ == '__main__':
    
    plt.close('all')
    plt.style.use('ggplot')
    
     
    xmax = 5
    ymax = 1
    tmax = 20
    scheme = 2
    p = 2
    jmaxes = [10,20,40,80]
    
    imaxes = [25,50,100,200]
    ind = 3

    imax = imaxes[ind]
    jmax = jmaxes[ind]

    fs = 20.
    fw = "bold"
    
    """
    cols = np.linspace(0,1,imax+1)
    rows = np.linspace(0,1,imax+1)
    
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
    ax1.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(imax),
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
    ax2.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(imax),
                  fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    
    """
    # ========================
    
    """
    xpos,ypos,csol = read_files(imax,jmax,xmax,ymax,0,scheme,p)
    
    
    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(csol),np.max(csol),20)
    CS = ax3.contour(xpos,ypos,csol, levels, colors='k')
    ax3.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(xpos,ypos,csol, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title('Contour Plot of the Solution Mesh for the Energy Equation, Meshsize = {}, t = 0'.format(imax),
                  fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\overline{T}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20) 
    """
    
    l2_norms = []
    for each_size in np.arange(0,4):
        
        imax = imaxes[each_size]
        jmax = jmaxes[each_size]
        
        xpos,ypos,csol = read_files(imax,jmax,xmax,ymax,tmax,scheme,p)
        
        Tw = analytic_sol(ypos)
        newTw = np.transpose(np.tile(Tw,(imax,1)))    
        
        error = newTw - csol
        
        
        l2_norms.append(calc_l2_norm(error,imax,jmax))
        
        if each_size == 0:
            fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
            levels = np.linspace(np.min(csol),np.max(csol),20)
            CS = ax4.contour(xpos,ypos,csol, levels, colors='k')
            ax4.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
            contour_filled = ax4.contourf(xpos,ypos,csol, levels)
            cbar = plt.colorbar(contour_filled)
            ax4.set_xlabel('x', fontsize=fs, fontweight=fw)
            ax4.set_ylabel('y', fontsize=fs, fontweight=fw)
            ax4.set_title('Contour Plot of the IE Solution for the Energy Equation, {} X {}, t = {}'.format(imax,jmax,tmax),
                          fontsize=fs,fontweight=fw)
            ax4.tick_params(axis='both',which='major',labelsize=20)
            cbar.set_label(r'$\overline{T}$',size=fs+5)
            cbar.ax.tick_params(labelsize=20) 
            
            fig_name = 'mech510_p3_ie_computed.png'
            plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches="tight", pad_inches=0.1,
                frameon=None)
                
                
            fig5,ax5 = plt.subplots(1,1,figsize=(16,12))
            levels = np.linspace(np.min(newTw),np.max(newTw),20)
            CS = ax5.contour(xpos,ypos,newTw, levels, colors='k')
            ax5.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
            contour_filled = ax5.contourf(xpos,ypos,newTw, levels)
            cbar = plt.colorbar(contour_filled)
            ax5.set_xlabel('x', fontsize=fs, fontweight=fw)
            ax5.set_ylabel('y', fontsize=fs, fontweight=fw)
            ax5.set_title('Contour Plot of the Steady-State Solution for the Energy Equation, {} X {}'.format(imax,jmax),
                          fontsize=fs,fontweight=fw)
            ax5.tick_params(axis='both',which='major',labelsize=20)
            cbar.set_label(r'$\overline{T}$',size=fs+5)
            cbar.ax.tick_params(labelsize=20) 
            
            fig_name = 'mech510_p3_ie_ana.png'
            plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches="tight", pad_inches=0.1,
                frameon=None)
       
    
    fig6,ax6 = plt.subplots(1,1,figsize=(16,12))
    
    l2_norms.sort()
    l2_norms = l2_norms[::-1]
    
    log_mesh = np.log10(imaxes)
    log_norms = np.log10(l2_norms)
    
    reg = np.polyfit(log_mesh, log_norms, 1)
    fit = np.polyval(reg, log_mesh)
    
    ax6.plot(log_mesh,log_norms,'or')
    ax6.plot(log_mesh,fit,'-b')
    
    print(reg)
    
    ax6.set_xlabel(r'log$_{10}$(Meshsize)', fontsize=fs, fontweight=fw)
    ax6.set_ylabel(r'log$_{10}$(L$_2$-Norm)', fontsize=fs, fontweight=fw)
    ax6.set_title(r'Log-Log Plot of IE L$_2$-Norms for Different Mesh Sizes', fontsize=fs, fontweight=fw)
    ax6.tick_params(axis='both',which='major',labelsize=20)
    ax6.set_xlim([0.8,3.0])
    fig_name = 'mech510_p3_ie_l2.png'
    
    plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
   
   
    
   
   
   
   
  
  

    