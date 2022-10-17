# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:58:37 2020

@author: ssm2617
"""

from oslo_model import * # PICKLE DOESNT LIKE THIS 
import matplotlib.pyplot as plt 
import pickle 
import scipy.optimize as op


def FindSteady(sizes, no_reps, plot = True):  
    '''
    > Pass an array of simulation sizes (parameter 'sizes') to find the number 
    of additions until the simulation has reached steady state (cutoff time). 
    > Each simulation is run a number (parameter 'repetitions') of times 
    and the time to steady state averaged over these repetitions.
    > The range of the cutoff time for each simulation size is also recorded 
    for the plot but only the average cutoff time is returned.
    
    Parameters
    ----------
    > sizes : list or array-like, 1d, dtype=int 
        list or 1d-array of system sizees to run the simulations for 
    > repetitions : int
        numer of times the simulation is run for each system size 
        
    Returns
    ----------
    > ave : list, dtype=int 
        list where each element is the number of avalanches of eacch system 
        size, averaged over (repetitions) number of runs 
    '''
    times = []    
    ave = []
    for i in sizes:
        rep_arr = [] # will contain all cutoff times for a single system size 
        for j in trange(no_reps, desc=f'L={i}'):
            s = OsloSim(i)
            s.RunTilSteady()
            rep_arr.append(s.steady_time)
        times.append(rep_arr) # recorded for range later -- not necessary to store here but could be useful if function is modified later    
        ave.append(np.average(rep_arr)) # average cutoff tiem for simulation size 
        
        
    if plot:
        min_err,max_err = [],[]
        for i in range(len(sizes)):
            min_err.append(ave[i] - np.min(times[i])) # for errorbars 
            max_err.append(np.max(times[i]) - ave[i]) 
        plt.figure()
        plt.errorbar(sizes,ave,yerr=[min_err,max_err], fmt = 'x')
        plt.xlabel('L')
        plt.ylabel('$<t_c>$')
        
    return ave
   
    
def SmootH(no_reps, size, run_length):
    '''
    > Initialise a number of simulations of same system size.
    > Run for a specified number of steps
    > Return total height of the pile over time, averaged over the 
    repeated simulations
    ------------
    PARAMETERS 
    > no_reps : int 
    Number of repeats to average over 
    > size : int 
    System size 
    > run_length : int 
    Number of grains to add to the system 
    ------------
    RETURNS 
    ave_height : list 
    '''
    heights = np.zeros([no_reps,run_length])
    
    for i in trange(no_reps):
        s = OsloSim(size)
        s.Run(run_length)
        heights[i] = s.heights  
    
    ave_height = [np.average(heights[:,j]) for j in range(run_length)]
    return ave_height
   
    
def HeightColl(height_arr, L):
    '''
    > Takes in array of heights from a simulation and generates data collapse
    of height/system size as a function of time 
    -----------
    PARAMETERS
    > height_arr : list or 1d-array
    list of values corresponding to a pile's total height over time 
    > L : int 
    system size 
    -----------
    RETURNS 
    > x : list of floats
    time/L^2 -- the x variable that leads to data collapse 
    > y : list of floats
    height/L -- the y variable that leads to data collapse 
    '''
    t_max = len(height_arr)
    height_arr = np.array(height_arr)
    t = np.linspace(1,t_max, len(height_arr))
    y = height_arr / L
    x = t / (L**2)
    
    plt.plot(x,y, linewidth = 0.6, label = L)
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel('t/L^2')
    plt.ylabel('L^-1 . h(t;L)')
    plt.legend()
    return x,y   

def PrepSteadySim(sim_sizes=[4,8,16,32,64,128,256,512], index=0):
    '''
    > Initialises and runs simulations of passed size until they have reached steady state
    > Uses Pickle to write out to file 
    > Simulations returned in array 
    -----------
    PARAMETERS 
    > sim_sizes : list of integers 
    Contains the system sizes to run simulations for 
    > index : integer 
    Append this integer to file name in case user wants to repeat preparation 
    of simulations 
    -----------
    RETURNS 
    > sim_arr : list of OsloSim objects 
    Contains simulations run until steady state 
    '''
    sim_arr = []
    for i in tqdm(sim_sizes):
        s = OsloSim(i)
        s.RunTilSteady()
        s.Run(100000)
        sim_arr.append(s)
        f = open('./pickles/steady_sims/all_steady_sims_%s.pkl'%index, 'wb')
        with f as output:
            pickle.dump(sim_arr,output,pickle.HIGHEST_PROTOCOL)
        f.close()
    return sim_arr

def PlotDiff(l_arr, h_arr):
    '''
    > h_arr[i] is steady-state height of simulation of size l_arr[i]
    > Reveals corrections to scaling (DIFFerence to naive) by dividing h by L 
    -----------
    PARAMETERS 
    > l_arr : list of integers 
    System size 
    > h_arr : list of floats 
    total pile height after reaching steady state 
    '''
    l_arr = np.array(l_arr)
    h_arr = np.array(h_arr)
    x = l_arr
    y = h_arr/l_arr
    plt.plot(x,y, 'x')
    plt.xlabel('L')
    plt.ylabel('$L^{-1} <h(t;L)>_t$')
    
    
def ScalH(l,a0,a1,w1):
    return a0*l*(1-a1*l**(-w1))
    
def ScalDiff(l,a0,a1,w1)  :
    return a0*(1-a1*l**(-w1))


#%%------------------------------- Question 2g ------------------------------

def SteadyH_L(sims, run_length):
    '''
    > Passed simulations and how many timesteps to consider
    > Finds and returns first and second moment of distribution 
    -----------
    PARAMETERS 
    > sims : list or 1d-array 
    List of OsloSim objects run for some time after reaching steady state 
    > run_length : int 
    The time for which the simulations were run after reaching steady-state 
    -----------
    RETURNS 
    > ave_height : list of floats
    List of average height of each simulation 
    > ave_height_sq : list of floats 
    List of average height^2 of each simulation 
    '''
    ave_height = np.array([])
    ave_height_sq = np.array([])
    for i in tqdm(sims):
        ave_height = np.append(ave_height, np.average(i.heights[-run_length:]))
        ave_height_sq = np.append(ave_height_sq, np.average((i.heights[-run_length:])**2))
    return ave_height, ave_height_sq

def H_SigH(sims,run_length):
    '''
    > Calculates sigma of steady-state height distribution from first and 
    second moments, given as sigma^2 = sqrt(<h^2> - <h>^2)
    -----------
    PARAMETERS 
    > sims : list or 1d-array 
    List of OsloSim objects run for some time after reaching steady state 
    > run_length : int 
    The time for which the simulations were run after reaching steady-state
    -----------
    RETURNS 
    > h : list of floats
    List of average height of each simulation 
    > sigma_sq : 1d-array 
    Array of sigma^2 for each simulation 
    '''
    h,h_sq = SteadyH_L(sims,run_length)
    return h, np.sqrt(h_sq - h**2)
    

def HDist(sim,no):
    '''
    > Take sim and range (from end) to find height distribution (only consider after tc)
    > Iterate through heights array and tick up relevant heights 
    > Returns array of distribution 
    -----------
    PARAMETERS 
    > sim : OsloSim object 
    the simulation to calculate the distribution for 
    > no : int 
    The number of positions from the end to calculate the ditribution for 
    -----------
    RETURNS 
    > hdist : 1d-array 
    array of height distribution 
    '''
    heights = sim.heights[-no:]
    # heights = sim.avalanches[-no:]
    max_h = int(max(heights))+1 # +1 because not counting 0 
    h_dist = np.zeros(max_h) # position denotes pile height 
    for i in range(no): # iterate through distribution array 
        h_dist[int(heights[i])] +=1       
    return h_dist/no


def AllHDist(sim_l, sims, time, plot = True ):
    '''
    > Takes array of system size, simulations themselves, and what timelength to 
    analyse over. 
    > Plots distribution for each system size and fits Gaussian to it 
    > Writes out all distributions as well as fit parameters 
    -----------
    PARAMETERS 
    > sim_l: list or 1d-array 
    List of system sizes for simulations 
    > sims : OsloSim object 
    Simulations to get distributions for 
    > time : int 
    Number of heights to sample 
    > plot : Boolean 
    If True, plot distributions. If False, do not plot distributions 
    -----------
    RETURNS 
    > all_h_dist : list 
    List contatining height distribution of each simulation 
    > Gaussian_params : list 
    List where each element corresponds to a simulation and contains bset-fit 
    parameters of that simulation's height distribution to a Gaussian.
    The paratmers are as follows: x0 (mean) and sigma (spread)
    '''
    all_h_dist = []
    params = []
    Gaussian_params = []
    for i in range(len(sims)):
        h_dist = HDist(sims[i],time)
        all_h_dist.append(h_dist)
 
        f = open('./pickles/h_data/h_dist.pkl','wb') # write out every time a system sieze has been analysed
        with f as output:
            pickle.dump(all_h_dist,output,pickle.HIGHEST_PROTOCOL)
        f.close()        
        x = np.arange(0,len(h_dist))
        y = h_dist
        params = op.curve_fit(Gauss,x,y,p0=[1.7*sim_l[i],1]) # from a0 and width is order of unity
        
        start = 0
        for j in range(len(h_dist)): # cutting out the non-zero heights
            if start == 0 and h_dist[j] != 0:
                start = j - 2 # have a couple of zeros 
        end = len(h_dist)
        x = np.arange(start,end)
        y = h_dist[start:end]
        Gaussian_params.append(params)
        g = open('./pickles/h_data/h_params.pkl','wb') # dump the gaussian h params 
        with g as output:
            pickle.dump(Gaussian_params,output,pickle.HIGHEST_PROTOCOL)
        g.close()
        
        if plot:
            plt.plot(x,y,'x', color = 'c', markersize = 6)
            x_fit = np.linspace(start,end,10000)    
            y_fit = [Gauss(k,params[0][0],params[0][1]) for k in x_fit]
    

            plt.plot(x_fit,y_fit, label = 'L = %s'%sim_l[i])
            plt.legend()
            plt.xscale('log')
            plt.grid()
            plt.xlabel('$h$')
            plt.ylabel('$P(h;L)$')
        
    return all_h_dist , Gaussian_params

def ExtractParams(all_params):
    '''
    Extracts h_mean and h_sig from all_params array 
    '''
    no_sims = len(all_params)
    h_mean = [all_params[k][0][0] for k in range(no_sims)]
    h_sig = [all_params[k][0][1] for k in range(no_sims)]
    return h_mean,h_sig

def CheckNormalised(all_d_hist):
    return [sum([i for i in all_d_hist[j]]) for j in range(len(all_d_hist))]

def HDistCollapse(h_dists,h_mean,h_sig):
    '''
    > Plot data collapse of pile hight as a Gaussian
    '''
    plt.figure()
    all_x = np.array([])
    all_y = np.array([])
    for i in trange(len(h_dists)): # i labels simulation length 
        h = np.arange(0,len(h_dists[i])) 
        for j in range(len(h_dists[i])): # for each simulation 
            x = []
            y = []            
            x.append( (h[j] - h_mean[i])/h_sig[i] )
            y.append(h_sig[i] * h_dists[i][j])            
            plt.plot(x,y,'.')
            all_x = np.append(all_x,x)
            all_y = np.append(all_y,y)
    x_plot = np.linspace(-4,4,10000)
    y_plot = [Gauss(i,0,1) for i in x_plot]
    plt.plot(x_plot,y_plot,'-', linewidth = 0.8, label = 'G(x; 0,1)')
    plt.legend() 
    
    
    x0,sig = op.curve_fit(Gauss,all_x,all_y)[0]
    y_fit = [Gauss(i,x0,sig) for i in x_plot]
    plt.plot(x_plot,y_fit,'-', linewidth = 0.8, label = 'G(x;%.4G,%.4G)'%(x0,sig))
    plt.legend()
    plt.xlabel('\frac{h - <h>}{\sigma}')
    plt.ylabel('\sqrt{\sigma} P(h;L)')

def Gauss(x,x0,sig):
    return 1/(np.sqrt(2*np.pi)*sig) * np.exp ( - ((x-x0)**2) /(2*(sig)**2) )



def DistributionQuestions():
    '''
    > Example of what order to use functions in order to plot data collapse and
    height distributions for data from first run of simulations  
    ----------
    RETURNS 
    > all_h_dist : list of floats 
    List of height distributions for each simulation 
    > Gaussian_params : list of floats 
    List of best-fit parameters for each simulation fitted to a Gaussian 
    
    '''
    sims,sims_1,sims_2, sims_3, sims_4 = GetSims()
    l_arr = [2**n for n in range(2,10)]
    h_arr = [sims[i].heights[-100000:] for i in range(8)] # Data Collapse for height over time 
    [HeightColl(h_arr[i],l_arr[i]) for i in range(2,10)]
    PlotDiff(l_arr,h_arr)
    all_h_dist , all_params = AllHDist(l_arr,sims,100000, plot = True )
    h_mean,h_sig = ExtractParams(all_params)
    HDistCollapse(all_h_dist,h_mean,h_sig)
    return all_h_dist, [h_mean,h_sig] 