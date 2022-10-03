# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 11:59:28 2019

@author: ssm2617

Passed parameters required to initialise Sim and Ball class
Initialises Simulation with passed initial conditions
Runs for 2000 frames to thermalise
Measure pressure

Increase number of balls in Container by 15 each simulation
Run for 15 data points
"""

from thermoCalc import *
import numpy as np 
import matplotlib.pyplot as plt 
from tqdm import trange

def pN(cont,n,m,r,v_bottom,v_top):
    '''
    See module docstring for information on purpose
    
    Parameters:
        cont = Container radius
        n = number of balls
        m - mass of balls
        r = radius of balls
        v_bottom = lower bound of velocity component speed
        v_top = upper bound of velocity componenet speed
    '''
    assert isinstance(n,int)
    assert isinstance(cont,float)    
    assert isinstance(m,float)
    assert isinstance(r,float)
    s=Sim(cont)
    pres_array=[]
    n_array=[]
    
    for i in trange(0,15):
        n_var = n+ 15*i
        s.generateBalls(n_var, m, r, v_bottom, v_top)
        s.Run(2000)
        n_array.append(n_var)        
        p=Pres(s)
        pres_array.append(p)
    return n_array, pres_array
    