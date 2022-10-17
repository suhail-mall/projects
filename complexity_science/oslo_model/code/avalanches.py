# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:10:03 2020

@author: suhai
"""

from oslo_model import * 
from logbin230119 import *
import scipy.optimize as op 
from heights import * 
import scipy as sp 


def CombineAvaArr():
    '''
    Takes avalanche data from pre-run sims Each sims_i contains data for each system length. 
    Returns array where each element is a list of avalanches
    '''
    sims,sims_1,sims_2, sims_3, sims_4 = GetSims()
    ava_arr = []
    for i in trange(len(sims)):
        arr = np.append(sims[i].avalanches[-100000:], sims_1[i].avalanches[-100000:])
        arr = np.append(arr,  sims_2[i].avalanches[-100000:])
        arr = np.append(arr,  sims_3[i].avalanches[-100000:])
        arr = np.append(arr,  sims_4[i].avalanches[-100000:])
        arr_int = [int(j) for j in arr]
        ava_arr.append(arr_int)
    return ava_arr    


def AllLogBinPlot(ava_arrs, l_arr, scale=1.5, zeros=False):
    '''
    > Input array of avalanche data for one sim
    > Returns binned distribution for all lengths (in list)
    ----------
    PARAMETERS 
    > ava_arrs : list where each element contains avalanche arrays 
    > l_arr : list of floats 
    Array containing system size for each simulation
    > scale : float, >1
    Parameter for logbinning 
    > zeros : Boolean 
    True if require distributionm False if just for plotting 
    ----------
    RETURNS 
    > all_x : list where each element contains logbin x array for a simulation
    > all_y : list where each element contains logbin y array for a simulation
    '''
    all_x = []
    all_y = []
    
    
    for i in trange(len(ava_arrs)):
        arr = ava_arrs[i]
        x,y = logbin(arr,scale, zeros)
        all_x.append(x)
        all_y.append(y)
        
        plt.plot(x,y,label = 'L = %s'%l_arr[i], linewidth = 0.7)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()    
    plt.grid()
    plt.xlabel('s')
    plt.ylabel('P(s;L)')
    return all_x,all_y 


def AvaDistCollapse(l_arr, all_x, all_y, tau_s = 1.55,D = 2.25):
    '''
    > Calculates parameters for data collapse 
    > Returns x and y for data collapse     
    > x = s/L^D , where s is avalanche size 
    > y = s^tau_s . P(s;L)
    ----------
    PARAMETERS 
    > l_arr : list of floats 
    List of simulation system sizes 
    > all_x : list where each element contains logbin x array for a simulation
    > all_y : list where each element contains logbin y array for a simulation
    > tau_s : float 
    Data Collapse parameter 
    D : float 
    Data Collapse parameter 
    ----------
    RETURNS 
    > collapse_x : 1d-array 
    x values for data collapse
    > collapse_y : 1d-array 
    y values for data collapse 
    '''
    plt.figure()
    for i in range(len(all_x)):
        collapse_x = np.array(all_x[i]) / (l_arr[i])**D # turn into array so can divide directly 
        collapse_y = np.array(all_y[i]) * np.array(all_x[i])**(tau_s)
        plt.plot(collapse_x,collapse_y, label = 'L = %s'%l_arr[i], linewidth = 0.7)
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    plt.xlabel('s/ L^D')
    plt.ylabel('s^tau_s . P(s;L)')
    return collapse_x, collapse_y


def fColl(s,a, tau):
    '''
    > Data Collapse function 
    > s is avalanche size 
    > a is multiplicative factor 
    > tau is power law 
    '''
    return a*s**-tau

def FitCollapse(all_x,all_y, start,end):
    '''
    Trial different regions of the plot to fit for tau_s
    '''
    x = all_x[-1][start:end]
    y = all_y[-1][start:end]
    a, tau_s= op.curve_fit(fColl,x,y,p0=[1.55])[0]
    x_plot = np.linspace(1,1e5,int(1e5))
    y_plot = [fColl(i,tau_s) for i in x_plot]
    plt.plot(x_plot,y_plot, '--',label = '%s'%tau_s)
    plt.legend()
    return a, tau_s    

def Moments(ava_sim, k_max):
    s_k = np.zeros(k_max) 
    T = len(ava_sim)
    for i in ava_sim:
        for j in range(len(s_k)):
            s_k[j] += i**(j+1)
    s_k/=T
    return s_k


def AllMoments(ava_arr,k_max):
    '''
    Calculates all moments
    Returned as arrays of moments for each L
    '''
    moments = []
    no_sims = len(ava_arr)
    for i in trange(no_sims):
        moments.append(Moments(ava_arr[i],k_max))
    return moments


def MomentPlot(l_arr,all_mom,k_max): 
    '''
    Takes all L to plot, all moments, and maximum moment
    Returns k_max lots of moments as function of L
    '''
    lines=[]
    for k in trange(k_max):
        x=[]
        y = []
        for l in range(len(all_mom)):
            x.append(l_arr[l])
            y.append(all_mom[l][k])
        plt.plot(x,y,'x',color = 'c')
        lines.append(y)
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    return lines

def MomFunc(L,a,b):
    return a*(L**b) 


def MomFuncLog(L,a,b):
     return a + L*b
 
def GradLinesLog(l_arr, lines):
     '''
     Find Grad of <s^k> against L lines 
     '''
     grads = []
     x = np.linspace(1,512,512)
     l_log = np.log(np.array(l_arr))
     lines_log = np.log(np.array(lines))
     for i in range(len(lines)):
         a,b = sp.polyfit(l_log,lines_log[i],1)
         grads.append(a)
         plt.plot(x, np.exp(b) *x**a, label = 'k = %s'%(i+1))
     plt.legend()
     return grads

def ScalingExpo(grads):
    '''
    Given gradient of moment plot, finds scaling exponents
    '''
    k_max = len(grads)
    k_arr = np.arange(1,k_max+1)
    (D,intercept), cov = sp.polyfit(k_arr,grads,1, cov = True)
    tau_s = 1 - intercept/D
    
    plt.plot(k_arr,grads, 'x')
    plt.grid()
    plt.xlabel('k')
    plt.ylabel('b = D(1 + k - tau_s)')
    x = np.linspace(0,12,100)
    y = [(D*i + intercept) for i in x]
    plt.plot(x,y,'-')
    return D, tau_s
    
def corr(x,a,b,c,d):
    return  b + a*x + c*x**-d   

def AvalancheExample():
    '''
    Example of what order the functions are run to find values of D and tau_s
    '''
    l_arr = [2**n for n in range(2,10)]
    ava_arr  = CombineAvaArr()
    all_x, all_y = AllLogBinPlot(ava_arr,l_arr)
    d_temp,tau_temp = FitCollapse(all_x,all_y,20,50)
    x_coll,y_coll = AvaDistCollapse(l_arr,all_x,all_y,tau_temp,d_temp) # plot with estimate values 
    
    all_x_0,all_y_0 = AllLogBinPlot(ava_arr,l_arr,zeros = True)# for moment analysis 
    mom = AllMoments(ava_arr,12)
    lines = MomentPlot(l_arr,mom,12) # plot moments (L)
    grads = GradLinesLog(l_arr,lines) # find gradients of moment plots 
    D,tau_s = ScalingExpo(grads)
    return D,tau_s