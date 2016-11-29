#!/usr/bin/env Python3.5

import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol
from sympy import sin,sqrt,Eq
from sympy.plotting import plot_implicit

def read_files(meshsize, xmax, tmax, cfl, ic):
    fname = "compute_" + str(meshsize) + "_xmax_" + str(xmax) + "_tmax_" + str(tmax) + "_cfl_" + str(cfl) + "_ic_" + str(ic) + ".txt"
    print(fname)
    csolution = []
    with open(fname, 'r') as cfile:
        for line in cfile:
            csolution.append(float(line.split()[0]))
            
    csolution = np.array(csolution)

    dx = float(xmax)/float(meshsize)
    xpos = np.array([dx*(i+0.5) for i in range(0,meshsize)])
    
    return(xpos,csolution)
    
    

def calc_l2_norm(error, meshsize):
    error_sum = np.sum(error**2)
    l2_norm = np.sqrt(error_sum/meshsize)
    
    return(l2_norm)
    
def calc_contour(aval,bval):
    sigma = np.sqrt((1+aval+(aval**2-bval**2)/2)**2 + (bval + aval*bval)**2)
    
    if (np.abs(sigma-1) < 1e-3):
        return(True, np.array([aval,bval]))

    else:
        return(False, [])

if __name__ == '__main__':
     
    xmax = 20
    tmax = 8
    cfl_list = [0.495,0.500,0.505,0.510][::-1]
    cfls = [int(np.round(each_cfl*1000)) for each_cfl in cfl_list]
    ic = 2
    
    meshsize = 500

    fs = 20.
    fw = "bold"
   
    plt.close('all')
        
    plt.style.use('ggplot')
    
    fig0,ax0 = plt.subplots(1,1,figsize=(16,12))
    
    l2_norms = []
    color_ind = 0
    colors = ['k','c','y','r'][::-1]
    for cfl in cfls:
        xpos,csolution=read_files(meshsize, xmax, tmax, cfl, ic)
        label_str = "CFL = " + "{:.3f}".format((cfl/1000.))
        color = colors[color_ind]
        color_ind = color_ind + 1
        
        if color == 'c':
            linewidth = 5
        else:
            linewidth = 1
        
        if color != 'k':
            ax0.plot(xpos, csolution, label=label_str, color=color, linewidth=linewidth)
        else:
            ax0.plot(xpos, csolution, label=label_str, color=color,linestyle='--',linewidth=linewidth)
    
    ax0.set_xlabel(r'$x_i$', fontsize=fs, fontweight=fw)
    ax0.set_ylabel(r'$\overline{T_i}$', fontsize=fs, fontweight=fw)
    ax0.set_title("Time Series of the Computed Solutions for Different CFL", fontsize=fs, fontweight=fw)    
    ax0.legend(loc='upper left',fontsize=fs)
    
    ax0.set_ylim([-1.5,1.5])
    
    fig_name = 'mech510_p2_3.png'
    fig0.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    a_array = np.linspace(-2,0,1000)
    b_array = np.linspace(-1.8,1.8,10000)
    
    #fig1,ax1 = plt.subplots(1,1,figsize=(16,12)) 
    """
    for each_a in a_array:
        for each_b in b_array:
            isappend,vals = calc_contour(each_a,each_b)
            if isappend:
                ax1.plot(vals[0],vals[1],'or')
    """
   
    x,y = Symbol('x'),Symbol('y') # Treat 'x' and 'y' as algebraic symbols
    p1= plot_implicit(Eq(sqrt((1+x+(x**2-y**2)/2)**2 + (y + x*y)**2), 1),label='hello',legend=True)
    fig1, ax1 = p1._backend.fig, p1._backend.ax
    fig1 = plt.gcf()
    fig1.set_size_inches(16, 12, forward=True)
    #fig1.figure(figsize=(16,12))    
    ax1.set_xlabel(r'Re($\lambda \Delta t$)', fontsize=fs, fontweight=fw)
    ax1.set_ylabel(r'Im($\lambda \Delta t$)', fontsize=fs, fontweight=fw, rotation=0, position=(3,0.95))
    ax1.set_title("Eigenvalues of the 2nd-Order Upwind Scheme, Against the RK3 Stability Boundary", fontsize=fs, fontweight=fw)    
    
    ax1.set_xlim([-2.1,1])
    ax1.set_ylim([-2,2])
    phis = np.linspace(0, 2*np.pi, 1000)
    
    
    cfl_list = np.array(cfl_list)
    delta_x = xmax/meshsize
    
    dts = cfl_list*delta_x/2
    print(dts)
    u = 2
    color_ind = 0
    for dt in dts:
        color = colors[color_ind]
        xlambda = u/(2.*delta_x)*(4*np.cos(phis)-np.cos(2*phis)-3)*dt
        ylambda = u/(2.*delta_x)*(np.sin(2*phis) - 4*np.sin(phis))*dt
        cfl = cfls[color_ind]
        color_ind = color_ind + 1
        dt_str = "{:.4f}".format(dt)
        cfl_str = "{:.3f}".format(cfl/1000)
        label_str = "CFL = " + cfl_str + r", $\Delta t = $" + str(dt_str)
        ax1.plot(xlambda,ylambda,color=color,linewidth=3,label=label_str)
    
    ax1.legend(loc='lower right',fontsize=fs)

    fig_name = 'mech510_p2_4.png'
    fig1.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    
    
    