#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.ndimage

import sys

import warnings
warnings.filterwarnings("ignore")


def savefigure(fig_name,fignum):
    print("Saving to {}".format(fig_name))
    fignum.savefig(fig_name, dpi=100, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches="tight", pad_inches=0.1,
        frameon=None)
    
    return(None)

def calc_stream(jarray,dx):
    sarray_out = np.zeros(np.shape(jarray))
    
    # first i index of streamfunction from ghost cell
    sarray_out[:,0] = jarray[:,0]*dx/2.
    
    for i in range(1,np.shape(jarray)[1]-1):
        sarray_out[:,i+1] = sarray_out[:,i] - jarray[:,i]*dx
        
    return(sarray_out)
    
def calc_stream1(uarray,dy):
    sarray = np.empty_like(uarray)
    
    # first i index of streamfunction from ghost cell
    sarray[0,:] = uarray[0,:]*dy/2.
    
    for j in range(1,np.shape(uarray)[0]-1):
        sarray[j+1,:] = sarray[j,:] + uarray[j,:]*dy
        
    return(sarray)

def calc_stream2(uarray,varray,dy,dx):
    sarray = np.empty_like(uarray)
    sarray1 = np.empty_like(uarray)
    sarray2 = np.empty_like(uarray)
    
    # first i index of streamfunction from ghost cell
    sarray1[0,:] = 0
    sarray1[:,0] = 0
    
    sarray2[0,:] = 0
    sarray2[:,0] = 0
    
    sarray1[:,-1] = 0
    sarray1[-1,:] = 0
        
    sarray2[:,-1] = 0
    sarray2[-1,:] = 0

    
    
    
    for j in range(1,np.shape(uarray)[0]-1):
        sarray1[j,1:-1] = np.sum(uarray[0:j+1,1:-1]*dy)
     
    for i in range(1,np.shape(varray)[1]-1):
        sarray2[1:-1,i] = np.sum(-varray[1:-1,0:i+1]*dx)
        
    sarray = sarray1 + sarray2   
    
    return(sarray)
    
if __name__ == '__main__':
    
    
    argus = sys.argv
    
    if len(sys.argv) < 10:
        print("meshsize,ymax,tlength,beta,apres,omega,Tw,dt,tol")
        exit()
    
    meshsize = int(sys.argv[1])
    ymax = int(sys.argv[2])
    tmax = int(sys.argv[3])
    beta = int(sys.argv[4])
    apres = int(sys.argv[5])
    omega = int(sys.argv[6])
    Tw = int(sys.argv[7])
    dt = float(sys.argv[8])
    tol = int(sys.argv[9])
        
    #plt.close('all')
    plt.style.use('ggplot')
    
    timings = []
    

    fs = 20.
    fw = "bold"

    ufilename = "umesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"
    vfilename = "vmesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"
    pfilename = "pmesh_cn"+str(meshsize)+"_xmax_1_ymax_"+str(ymax)+"_t_"+str(tmax)+"_p_3_"+str(beta)+"_"+str(apres)+"_"+str(omega)+"_"+str(Tw)+".txt"

    uarray = []
    varray = []
    parray = []
    
    with open(ufilename, 'r') as ufile:
        for line in ufile:
            num_list = [float(num) for num in line.split()]
            uarray.append(num_list)    
    with open(vfilename, 'r') as vfile:
        for line in vfile:
            num_list = [float(num) for num in line.split()]
            varray.append(num_list)
    with open(pfilename, 'r') as pfile:
        for line in pfile:
            num_list = [float(num) for num in line.split()]
            parray.append(num_list)

    uarray = np.flipud(np.array(uarray))
    varray = np.flipud(np.array(varray))
    speed = np.sqrt(uarray**2 + varray**2)
    
    parray = np.flipud(np.array(parray))
    
    filt = 1
    # Pressure and velocity plot
    speed = scipy.ndimage.zoom(speed,filt)
    uarray = scipy.ndimage.zoom(uarray,filt)
    varray = scipy.ndimage.zoom(varray,filt)
    
    imax = meshsize
    dx = dy = 1./imax
    
    jmax = (ymax/10)/dy
    
    cols = np.linspace(0,1,imax*filt+1)
    rows = np.linspace(0,ymax/10,jmax*filt+1)
    
    
    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = -(rows + np.diff(rows)[0]/2)[:-1]


    fig1,ax1 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(speed),np.max(speed),80)
    #CS = ax1.contour(cols_plot, rows_plot, parray, levels, colors='w',linestyles='dashed')
    #ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    #contour_filled = ax1.contourf(cols_plot, rows_plot, parray, levels)
    #CS = ax1.contour(cols_plot, rows_plot, speed, levels, colors='w',linestyles='dashed')
    #ax1.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    #contour_filled = ax1.contourf(cols_plot, rows_plot, speed, levels)
    #cbar = plt.colorbar(contour_filled)
    ax1.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax1.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax1.set_title(r'Steady-State Velocity Streamlines and Speeds, Tw = {}, dt = {}, Meshsize = {} X {}'.format(Tw,dt,imax,int(jmax)),fontsize=fs,fontweight=fw)
    ax1.tick_params(axis='both',which='major',labelsize=20)
    #cbar.set_label(r'$P_{computed}$',size=fs+5)
    #cbar.ax.tick_params(labelsize=20)
    
    cols = np.linspace(0,1,imax*filt+1)
    rows = np.linspace(0,ymax/10,jmax*filt+1)


    cols_plot = (cols + np.diff(cols)[0]/2)[:-1]
    rows_plot = -(rows + np.diff(rows)[0]/2)[:-1]
    
    #domain = 5*len(speed)/6
    domain = 0
    X, Y = np.meshgrid(cols_plot,rows_plot)
    lw = speed[domain:,:]/speed.max()
    
    print(np.shape(lw))
    plt.streamplot(X[domain:,:],Y[domain:,:],uarray[domain:,:],varray[domain:,:],color=lw,cmap=cm.jet,density=5,zorder=5)
    cbar = plt.colorbar()
    cbar.set_label(r'Speed',size=fs)
    cbar.ax.tick_params(labelsize=20)
    print(meshsize)
    #figname = "pcontour_velocities_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.png".format(pcontour_velocities_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.png)
    #savefigure(figname,fig1)
    
    dx = dy = 1./meshsize
    
    domain = 10
    restrict = 1*np.shape(varray)[1]/5
    
    
    sarray = calc_stream(varray,dx)
    sarray_plot = (np.abs(sarray[domain:,restrict:-restrict]))

    maxind = np.where(np.abs(sarray) == np.max(sarray_plot))
    print(maxind)

    fig2,ax2 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(sarray_plot),np.max(sarray_plot),80)

    contour_filled = ax2.contourf(cols_plot[restrict:-restrict], rows_plot[domain:], sarray_plot, levels)
    cbar = plt.colorbar(contour_filled)
    ax2.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax2.set_title(r'Steady-State Streamline Contours, Tw = {}, dt = {}, Meshsize = {} X {}'.format(Tw,dt,imax,int(jmax)),fontsize=fs,fontweight=fw)
    ax2.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\psi_{i,j}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    ax1.plot(cols_plot[maxind[1]],rows_plot[maxind[0]],'og',markersize=20)
 
    print(cols_plot[maxind[1]])
    print(rows_plot[maxind[0]])
    print(sarray[maxind[0],maxind[1]])
    print(np.max(sarray_plot))
    
    
    fname = "streams_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_1.txt".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize,domain)
    with open(fname,'w') as outfile:
        outfile.write(str(cols_plot[maxind[1]][0])+"\n")
        outfile.write(str(rows_plot[maxind[0]][0])+"\n")
        outfile.write(str(np.max(sarray_plot)))
        
    print("")
    
    domain = 2*len(varray)/5
    sarray1 = calc_stream(varray,dx)
    sarray_plot1 = (np.abs(sarray1[domain:,restrict:-restrict]))
    
    maxind = np.where(np.abs(sarray1) == np.max(sarray_plot1))
    print(maxind)

    fig3,ax3 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(sarray_plot1),np.max(sarray_plot1),80)
    #CS = ax2.contour(cols_plot, rows_plot[domain:], sarray, levels, colors='w',linestyles='dashed')
    #ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax3.contourf(cols_plot[restrict:-restrict], rows_plot[domain:], sarray_plot1, levels)
    cbar = plt.colorbar(contour_filled)
    ax3.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax3.set_title(r'Steady-State Streamline Contours, Tw = {}, dt = {}, Meshsize = {} X {}'.format(Tw,dt,imax,int(jmax)),fontsize=fs,fontweight=fw)
    ax3.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\psi_{i,j}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    ax1.plot(cols_plot[maxind[1]],rows_plot[maxind[0]],'og',markersize=20)
    
    print(cols_plot[maxind[1]])
    print(rows_plot[maxind[0]])
    print(sarray1[maxind[0],maxind[1]])
    print(np.max(sarray_plot1))
    
    fname = "streams_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_2.txt".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize,domain)
    with open(fname,'w') as outfile:
        outfile.write(str(cols_plot[maxind[1]][0])+"\n")
        outfile.write(str(rows_plot[maxind[0]][0])+"\n")
        outfile.write(str(np.max(sarray_plot1)))
    
    print("")
    
    domain = 7*len(varray)/8
    
    sarray2 = calc_stream(varray,dx)
    sarray_plot2 = (np.abs(sarray2[domain:,restrict:-restrict]))

    maxind = np.where(np.abs(sarray2) == np.max(sarray_plot2))
    print(maxind)
    
    fig4,ax4 = plt.subplots(1,1,figsize=(16,12))
    levels = np.linspace(np.min(sarray_plot2),np.max(sarray_plot2),80)
    #CS = ax2.contour(cols_plot, rows_plot[domain:], sarray, levels, colors='w',linestyles='dashed')
    #ax2.clabel(CS, colors='k', fmt='%6.4f', fontsize=16)
    contour_filled = ax4.contourf(cols_plot[restrict:-restrict], rows_plot[domain:], sarray_plot2, levels)
    cbar = plt.colorbar(contour_filled)
    ax4.set_xlabel('x', fontsize=fs, fontweight=fw)
    ax4.set_ylabel('y', fontsize=fs, fontweight=fw)
    ax4.set_title(r'Steady-State Streamline Contours, Tw = {}, dt = {}, Meshsize = {} X {}'.format(Tw,dt,imax,int(jmax)),fontsize=fs,fontweight=fw)
    ax4.tick_params(axis='both',which='major',labelsize=20)
    cbar.set_label(r'$\psi_{i,j}$',size=fs+5)
    cbar.ax.tick_params(labelsize=20)
    
    ax1.plot(cols_plot[maxind[1]],rows_plot[maxind[0]],'og',markersize=20)
    
    print(cols_plot[maxind[1]])
    print(rows_plot[maxind[0]])
    print(sarray2[maxind[0],maxind[1]])
    print(np.max(sarray_plot2))

    fname = "streams_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_3.txt".format(ymax,tmax,beta,apres,omega,Tw,dt,tol,meshsize,domain)
    with open(fname,'w') as outfile:
        outfile.write(str(cols_plot[maxind[1]][0])+"\n")
        outfile.write(str(rows_plot[maxind[0]][0])+"\n")
        outfile.write(str(np.max(sarray_plot2)))
        
    #plt.show()
    
   