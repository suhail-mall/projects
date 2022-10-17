# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 17:47:25 2022

@author: suhai
"""
from heights import * 
       
def GetHDistData():
    '''
    > Gets height data from pickled data files 
    > Also returns the simulations that correspond to those data 
    ----------
    RETURNS 
    > h_dists : list of floats
    list where each element is the height distribution for each simulation
    > h_params : list of floats
    list where each element is the Gaussian fit parameters for each simulation
    > sims : list of OsloSim objects 
    List of simulations that were run and saved, and who generated the 
    height data 
    '''
    f = open(r'../pickles/h_data/h_dist.pkl','rb')
    with f as input:
        h_dists = pickle.load(input)
    g = open('../pickles/h_data/h_params.pkl', 'rb')
    with g as input:
        h_params = pickle.load(input)
    h = open(r'../pickles/steady_sims/all_steady_sims.pkl','rb')
    with h as input:
        sims = pickle.load(input)   
    return h_dists, h_params , sims
    
def GetHAveData():
    '''
    > Retrieve average heights over time for L = 4 to L = 256
    ----------
    RETURNS 
    > h_ave : list of floats 
    list of average heights, each element corresponding to a different 
    system size 
    '''
    f = open('../pickles/steady_sims/ave_h_4to256.pkl','rb')
    with f as input:
        h_ave = pickle.load(input)
    f.close()
    return h_ave
    
def GetSims():
    '''
    Retrieve steady-state simulations for L = 4 to L = 512
    Have 5 repeats 
    Each sims_i contains simulations for each L 
    ---------
    RETURNS 
    sims : tuple 
    > Tuple conttaining 5 elements
    > Each element contains a list of 7 simulations of length 
    4, 8, 16, 32, 64, 128, 256
    '''
    sims = []
    for i in range(5):
        f = open(r'../pickles/steady_sims/all_steady_sims_%s.pkl'%i,'rb')
        with f as input:
            sims.append(pickle.load(input))
        f.close()
    
    return tuple(sims)

