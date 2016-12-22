#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.ndimage

import sys
import glob
import warnings
warnings.filterwarnings("ignore")


def calc_error(phi1,phi2,phi3,h1,h2,h3):
    
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
    

    print(p)
    print(ea21*100)
    print(eext21*100)
    print(GCI21*100)
    print(ea32*100)
    print(eext32*100)
    print(GCI32*100)
    print("")
    return(None)

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
    
    s180 =  ["streams_30_15277_100_0_100_100_0.2_14_180_10_1.txt","streams_30_15277_100_0_100_100_0.2_14_180_216.0_2.txt","streams_30_15277_100_0_100_100_0.2_14_180_472.5_3.txt"]
    s169 = ["streams_30_13567_100_0_100_100_0.2_14_169_10_1.txt","streams_30_13567_100_0_100_100_0.2_14_169_202.8_2.txt","streams_30_13567_100_0_100_100_0.2_14_169_443.625_3.txt"]
    s130 = ["streams_30_8300_100_0_100_100_0.2_14_130_10_1.txt","streams_30_8300_100_0_100_100_0.2_14_130_156.0_2.txt","streams_30_8300_100_0_100_100_0.2_14_130_341.25_3.txt"]
    s120 = ["streams_30_7152_100_0_100_100_0.2_14_120_10_1.txt","streams_30_7152_100_0_100_100_0.2_14_120_144.0_2.txt","streams_30_7152_100_0_100_100_0.2_14_120_315.0_3.txt"]
    s100 = ["streams_30_5101_100_0_100_100_0.2_14_100_10_1.txt","streams_30_5101_100_0_100_100_0.2_14_100_120.0_2.txt","streams_30_5101_100_0_100_100_0.2_14_100_262.5_3.txt"]
    s80 = ["streams_30_3379_100_0_100_100_0.2_14_80_10_1.txt","streams_30_3379_100_0_100_100_0.2_14_80_96.0_2.txt","streams_30_3379_100_0_100_100_0.2_14_80_210.0_3.txt"]
    s40 = ["streams_30_946_100_0_100_100_0.2_14_40_10_1.txt","streams_30_946_100_0_100_100_0.2_14_40_48.0_2.txt","streams_30_946_100_0_100_100_0.2_14_40_105.0_3.txt"]
    s20 = ["streams_30_265_100_0_100_100_0.2_14_20_10_1.txt","streams_30_265_100_0_100_100_0.2_14_20_24.0_2.txt","streams_30_265_100_0_100_100_0.2_14_20_52.5_3.txt"]
    meshsizes = [80,100,120,130,169,180]
    
    
    #sfiles1 = glob.glob('streams*1.txt')[::-1]
    #sfiles2 = glob.glob('streams*2.txt')[::-1]
    #sfiles3 = glob.glob('streams*3.txt')[::-1]
    
    sfiles1 = [s180[0],s120[0],s80[0]]
    sfiles2 = [s180[1],s120[1],s80[1]]
    sfiles3 = [s180[2],s120[2],s80[2]]
    
    dx1 = dy1 = 1/180.
    dx2 = dy2 = 1/120.
    dx3 = dy3 = 1/80.
    
    
    
    sx1 = []
    sy1 = []
    sa1 = []
    
    for sfile1 in sfiles1:
        print(sfile1)
        
        temp_arr = []
        with open(sfile1,'r') as f:
            for line in f:
                temp_arr.append(float(line))
                
        sx1.append(temp_arr[0])
        sy1.append(temp_arr[1])
        sa1.append(temp_arr[2])
        
    sx2 = []
    sy2 = []
    sa2 = []
    
    for sfile2 in sfiles2:
        print(sfile2)
        
        temp_arr = []
        with open(sfile2,'r') as f:
            for line in f:
                temp_arr.append(float(line))
                
        sx2.append(temp_arr[0])
        sy2.append(temp_arr[1])
        sa2.append(temp_arr[2])    
    
    
    sx3 = []
    sy3 = []
    sa3 = []
    
    
    for sfile3 in sfiles3:
        print(sfile3)
        
        temp_arr = []
        with open(sfile3,'r') as f:
            for line in f:
                temp_arr.append(float(line))
                
        sx3.append(temp_arr[0])
        sy3.append(temp_arr[1])
        sa3.append(temp_arr[2])
        

    

    h1 = np.sqrt(dx1*dy1)
    h2 = np.sqrt(dx2*dy2)
    h3 = np.sqrt(dx3*dy3)       
    
    print("")
    print("First vortex \n")
    calc_error(sx1[0],sx1[1],sx1[2],h1,h2,h3)
    calc_error(sy1[0],sy1[1],sy1[2],h1,h2,h3)
    calc_error(sa1[0],sa1[1],sa1[2],h1,h2,h3)
                
    print("")
            
    print("Second vortex \n")
    calc_error(sx2[0],sx2[1],sx2[2],h1,h2,h3)
    calc_error(sy2[0],sy2[1],sy2[2],h1,h2,h3)
    calc_error(sa2[0],sa2[1],sa2[2],h1,h2,h3)
                
    print("")
    
    print("Third vortex \n")
    calc_error(sx3[0],sx3[1],sx3[2],h1,h2,h3)
    calc_error(sy3[0],sy3[1],sy3[2],h1,h2,h3)
    calc_error(sa3[0],sa3[1],sa3[2],h1,h2,h3)
    
    