# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:28:28 2020

@author: suhailMall
"""

import numpy as np
from tqdm import *


class OsloSim:
    '''
    Initialise a simulation of the Oslo Model
    '''
    
    def __init__(self, length = 64, prob = 0.5):
        '''
        > Simulation class containing all relevant parameters and functions 
        > Can run a simulation of the Oslo Model using Run() and Trun() 
        functions 
        > Also store properties of the simulated sand pile, e.g. number of 
        avalanches (relaxations) per drive, height of pile after each drive,
        number of drives until pile is stable.
        > Properties are extracted and used in avalanches.py and heights.py 
        files 

        Parameters
        ----------
        length : int
             > Length of the system
             > default is 64
        prob : float, between 0 and 1
             > Thershold slope can be either 1 or 2 
             > probability of the threshold slope being set to 1
             > P(z_th=1)=p, P(z_th=2)=1-p
             > default is 0.5.
        '''
        self.L = int(length)
        self.p = prob    
       
        self.h       = np.zeros(self.L)
        self.z_th    = np.random.choice([1,2], p=[self.p, 1-self.p], size=self.L)
        self.max_end = 0 # The furthest distance a grain has tumbled -- don't need to check past this 
       
        self.avalanche_size = 0 # For records 
        self.no_drives      = 0
       
        self.steady      = False # Steady if all sites filled and grain tumbles out of system 
        self.steady_time = 0 # Number of drives to become 
        
        self.heights    = np.array([]) # Total height of pile after all sites are relaxed, each time you add a grain
        self.avalanches = np.array([]) # Total number of relaxations each time you add a grain 
        
    def Z(self, i):
        '''
        > slope at site i, given as h_{i+1}-h_i 
        > If i is the final site (i=L-1), return h_i (i.e. assume h=0 
        outside of system)
        
        Parameters
        ----------
        i : int 
            index representing site to check height at 
            
        Returns
        -------
        Z : int
            slope at site i
        '''
        if i == self.L - 1: # Final position 
            return self.h[i]
        return self.h[i] - self.h[i+1]
   
    def Drive(self):
        '''
        >Adds a grain at 0th site -- update height
        '''
        self.h[0] += 1        
       
    def RelaxSite(self, site):
        '''
        > Remove a grain from site i and add it to site (i+1)
        
        Parameters
        ----------
        site : int
            index representing position of site i 
        '''
        self.h[site] -= 1
        if site != self.L - 1: # only add move grain to right if not at end of system
            self.h[site+1] += 1
           
        new_th = np.random.choice([1,2],p = [self.p, 1 - self.p])
        self.z_th[site] = new_th
       
       
    def Stable(self, i):
        '''
        > checks if site i is stable or not 
        > If slope at site i is greater than threshold at site i, then 
        site is unstable 

        Parameters
        ----------
        i : int
            index representing which site to check 

        Returns
        -------
        stable : Boolean
            True if site is stable, False if not. 
        '''
        return not self.Z(i) > self.z_th[i] 
       
       
    def IterForward(self,start):
        '''
        > Given a starting point to check, relax this point and consider the 
        site to the right. 
        > Relax this point and continue checking the next site until 
        
        Parameters
        ----------
        start : int
            index representing site to start iterating forward from 
        '''
        i = start
        iterate_forward = True
         # counts forward until you hit a stable site. Resets everytime system is driven
        while iterate_forward:            
            if not self.Stable(i):
               
                self.RelaxSite(i) # relax the site and follow to the next
                self.avalanche_size += 1 # Count number of relaxations for this drive 
               
                if not self.steady and i == self.L -1: # Checking if steady and counting time to steady 
                    self.steady = True
                    self.steady_time = self.no_drives - 1
                   
                i+=1 # move to next site
                
                if i == self.L: # already relaxed final site now (next site is at end)
                    iterate_forward = False # at end of system  
            else:
                iterate_forward = False # stop iterating forward once you hit a stable site
 
        if i-1 > self.max_end: # after finished iterating forward, i is the final site relaxed 
            self.max_end = i-1  # check if this is furthest point  -- i-1 is final RELAXED point
     
           
    def IterBack(self, i):
        '''
        > Relax the site, which can cause the sites either side of it to become 
        unstable
        > Therefore, relax the site then iterate forward, "following the grain"
        until it becomes stable
        > Then consider the site to the left and do the same 
        > Repeat until
        Terminate when hit a stable point or start of system
        
        Parameters 
        ----------
        i : int 
            index representing site to start iterating backward from 
        '''
        iterate_back = True # Initialise True to start iterating backwards 
        while iterate_back: # While the starting site of an avalanche is not stable 
            if self.Stable(i): # if starting site is stable, stop iterating back 
                iterate_back = False                
            self.IterForward(i) # Relax the site and follow the grain forward 
            if i == 0: # If reach start of pile 
                iterate_back = False
            i -= 1 # Move backwards to previous site 
        return i + 1 # first stable point
   
    def UntilStable(self):
        '''
        Checks relevant sites until all are stable
        '''        
        j = 0 # start at 0
        while j <= self.max_end: # Until we hit the unaffected region that has seen no relaxations this drive  -- this is checked and updated every relaxation 
            j = self.IterBack(j) # iterate backwards, and then use this final stable point as the start for the next check 
            j+=1 # move to next site until all affected sites checked
        self.max_end = 0 # reset the max affected range after system is stable
        
        # Store total height of pile, number of avalanches etc for records 
        self.heights = np.append(self.heights,self.h[0])
        self.avalanches = np.append(self.avalanches,self.avalanche_size)
        self.avalanche_size = 0
        
    def Trun(self, no):
        '''
        > Run for specified time, with progress bar using tqdm  
        > Time given by number of grains added
        
        Parameters
        ----------
        no : int
            number of grains to add to system 
        '''
        
        for i in trange(no):
            self.Drive()
            if not self.steady: # If the system is not steady yet, increase count till steady. Stop counting when steady  
                self.no_drives += 1  
            self.UntilStable()
 
    def Run(self, no):
        '''
        > Run for specified time, without progress bar 
        > Time given by number of grains added
        
        Parameters
        ----------
        no : int
            number of grains to add to system 
        '''
        for i in trange(no):
            self.Drive()
            if not self.steady: # If the system is not steady yet, increase count till steady. Stop counting when steady 
                self.no_drives += 1            
            self.UntilStable()    
                   
    def RunTilSteady(self):
        '''
        > Runs until a grain leaves the system
        > This implies that all sites in the system are occupied and stable 
        
        Returns
        ----------
        steady_time : int
            number of drives until system has become steady 
        '''
        while not self.steady:
            self.Drive()
            self.no_drives += 1
            self.UntilStable()
        return self.steady_time
       
    def Check(self):
        '''
        > Array of booleans with each element (i) showing if the site (i)is 
        stable 
        
        Returns 
        ----------
        arr : array_like, 1 dimensional, length L 
            Array of booleans with each element (i) showing if the site (i)is 
            stable 
        '''
        arr = [self.Stable(i) for i in range(self.L -1)]
        return arr


