# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:49:26 2022

@author: suhai



Just build framework for tree stuff 
Add Bethe Lattice in some other file 
"""
import numpy as np 


class Tree:
    def __init__(self, p=0.5):
        self.p = p 
        
        self.parent_node = np.random.choice([0,1], p=[1-self.p,self.p])
        
        
        
        return 
    

class node:
    def __init__(self, value, left, right):
        self.value = self.some_func(value) # In case not just write value directly 
        self.left_node  = left #daughter nodes 
        self.right_node = right 
        
    def some_func(self, value):
        return value 