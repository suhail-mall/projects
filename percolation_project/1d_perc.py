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
        
        # NEED TO IMPLMN PERIODIC BOUNDARY CONDNS  !!!!!!! Affects first and last cluster 
        clust_dist = np.zeros(self.L+1) # Largest possible cluster size is L, start counting at 0 to match indexing
        clust_l    = 0 # cluster length 
        
        for site in self.sites:
            if site:
                clust_l += 1
            elif clust_l: # Site not active and current length!=0
                clust_dist[clust_l] += 1
                clust_l = 0
        clust_dist[clust_l] += 1 # If cluster contains final site
        return clust_dist / np.sum(clust_dist)
    
    
    
    
    
    
def theo_nsp(s, p, L):
    
    if 0<=p<=1:
        raise Exception ('Probability not in range [0,1]')
    if 1<=s<=L:
        raise Exception('Cluster size not in range [1,L]')
    s = int(s)
    L = int(L)
    if s == L:
        return p **L /L
    if s == L - 1: 
        return (1-p) * p **(s-1)
    else: 
        return (1-p)**2 * p**s 



