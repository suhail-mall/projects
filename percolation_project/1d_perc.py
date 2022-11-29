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
        clust_l    = 0 # cluster length 
        
        for site in self.sites:
            if site:
                clust_l += 1
            elif clust_l: # Site not active and current length!=0
                clust_dist[clust_l] += 1
                clust_l = 0
        clust_dist[clust_l] += 1 # If cluster contains final site
        return clust_dist / np.sum(clust_dist)
    
    
    
    