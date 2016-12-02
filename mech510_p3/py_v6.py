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
    
     
    xmax = 15
    ymax = 1
    tmax = 40
    scheme = 2
    p = 3
    jmaxes = [10,20,40,80]
    
    imaxes = [25,50,100,200]
    ind = 1

    imax = imaxes[ind]
    jmax = jmaxes[ind]

    fs = 20.
    fw = "bold"
    
   

    l2_norms = []
    
    max_bots = []
    grids = []
    for each_size in np.arange(0,4):
        
        imax = imaxes[each_size]
        jmax = jmaxes[each_size]
        
        grids.append(imax*jmax)
        
        
        dx = xmax/imax
        dy = ymax/jmax
        
        xpos,ypos,csol = read_files(imax,jmax,xmax,ymax,tmax,scheme,p)
        
        max_bot = (np.max(csol[0,:]))
        max_loc = np.where(csol[0,:] == max_bot)
        max_x = (max_loc[0]+0.5)*dx
        
        print(max_bot)
        print(max_loc[0])
        print(max_x)
        
        max_bots.append(np.max(max_bot/dy))
        
        
        if each_size == 3:
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
            
        
            fig_name = 'mech510_p3_ie_computed_7.png'
            plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches="tight", pad_inches=0.1,
                frameon=None)
            
            fig5,ax5 = plt.subplots(1,1,figsize=(16,12))
            
            gradients = csol[0,:]/dy
            
            ax5.plot(xpos,gradients,'-r')
            ax5.set_xlabel('x', fontsize=fs, fontweight=fw)
            ax5.set_ylabel('dT/dy', fontsize=fs, fontweight=fw)
            ax5.set_title('Temperature Gradient at the Bottom Wall, {} X {}, t = {}'.format(imax,jmax,tmax),
                          fontsize=fs,fontweight=fw)
            ax5.tick_params(axis='both',which='major',labelsize=20)
            
            fig_name = 'mech510_p3_ie_gradient_7.png'
            plt.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches="tight", pad_inches=0.1,
                frameon=None)
                       
                       
    max_bots = max_bots[1:]
    grids = grids[1:]
    
    
    phi1 = max_bots[0]
    phi2 = max_bots[1]
    phi3 = max_bots[2]
    
    h1 = grids[0]
    h2 = grids[1]
    h3 = grids[2]
    
    r21 = h2/h1
    r32 = h3/h2
    
    eps32 = phi3 - phi2
    eps21 = phi2 - phi1
    
    p = 1/np.log(r21)*np.abs(np.log(np.abs(eps32/eps21)))
    
    phiext21 = ((r21**p)*phi1 - phi2)/(r21**p-1)
    phiext32 = ((r32**p)*phi2 - phi3)/(r32**p-1)
    
    ea21 = np.abs((phi1 - phi2)/phi1)
    eext21 = np.abs((phiext21 - phi1)/phiext21)
    GCI21 = (1.25*ea21/(r21**p-1))
    
    ea32 = np.abs((phi2 - phi3)/phi2)
    eext32 = np.abs((phiext32 - phi2)/phiext32)
    GCI32 = (1.25*ea32/(r32**p-1))
    
    print("")
    print(p)
    print("")
    
    print(ea21*100)
    print(eext21*100)
    print(GCI21*100)
    
    print("")
    
    print(ea32*100)
    print(eext32*100)
    print(GCI32*100)
    
    
    
    
    """
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
    """
   
   
    
   
   
   
   
  
  

    