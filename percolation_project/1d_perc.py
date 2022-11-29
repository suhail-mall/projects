# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:02:55 2022

@author: suhai
"""
import numpy as np


class Sim:
    def __init__(self, length=64, p=0.5):
        self.L = length
        self.p = p 
        self.populate_sites()
        
    def populate_sites(self):
        self.sites = np.random.choice([0,1], p=[1-self.p, self.p], size=self.L)
        
    def cluster_dist(self):
        '''
        > Calculate the distribution n(s,p) by counting the number of clusters 
        and grouping them by their size 
        > Checked by removing normalisation and summing i*dist[i] over dist,
        and comparing to sum(Sim.sites)

        Returns
        -------
        clust_dist : array_like, 1_dimenstional, length Sim.L+1
            The normalised cluster distribution 
        '''
        
        clust_dist = np.zeros(self.L+1) # Largest possible cluster size is L, start counting at 0 to match indexing
        
        clust_l  = 0 # cluster length 
        first_cl = 0
        if self.sites[0]: # Store length of first cluster if start at pos 0 (periodic bound condns)
            for site in self.sites:
                if site:
                    first_cl+=1
                else:
                    break
        
        for site in self.sites[first_cl+1:]:
            if site:
                clust_l += 1
            elif clust_l: # Site not active and current length!=0
                clust_dist[clust_l] += 1
                clust_l = 0
        clust_dist[clust_l + first_cl] += 1 # If cluster contains final site, and add first cl due to bound condns
        return clust_dist / np.sum(clust_dist)
    
    
    
 #%% Theoretical Predcitions and Plots 
    
    
def nspl(s, p, L):
    '''
    > Calculated probability of a cluster of size s in a system of length L and 
    activation probability p
    > Denoted n(s,p,L)
    Parameters
    ----------
    s : int, [1,L]
        Size of cluster
    p : float, [0,1]
        Probability of each site being activated 
    L : int
        Length of system
    Returns
    -------
    n : float
        n(s,p,L)
    '''
    if not 0<=p<=1:
        raise Exception ('Probability not in range [0,1]')
    if not 1<=s<=L:
        raise Exception('Cluster size not in range [1,L]')
    s = int(s)
    L = int(L)
    if s == L:
        return p **L /L
    if s == L - 1: 
        return (1-p) * p **(s-1)
    else: 
        return (1-p)**2 * p**s 


#%% Testing 

def test_dist(L = 2**10, reps = 2000):
    for i in range(reps): # Do by removing normalisation in final cluster_dist line 
        s = Sim(L)
        d = s.cluster_dist()
        # print(i, sum(s.sites), sum([j*d[j] for j in range(len(d))]))
        print(i,  sum(s.sites) == sum([j*d[j] for j in range(len(d))]))
    
    

