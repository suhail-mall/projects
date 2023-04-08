# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:49:26 2022

@author: suhai



Just build framework for tree stuff 
Add Bethe Lattice in some other file 
"""
import numpy as np 
    
#%%

class Tree:
    def __init__(self, value, gen_p=0.5, vals=[0.5,0.5], val_p=[0.5,0.5]):
        self.value   = value
        self.subtree = self.gen_next(gen_p, vals, val_p)
        
    def gen_next(self, gen_p, vals, val_p):
        if_gen = np.random.choice(a=[True,False], p=[gen_p, 1-gen_p], size=2) 
        # Can replace this with a function to choose if_gen
        subtree = [Tree(self.value) for _ in range(np.sum(if_gen))] 
        # Operate on subtree here 
        return subtree 