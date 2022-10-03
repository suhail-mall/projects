# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:02:34 2019

@author: ssm2617

Passed parameters required to initialise Sim and Ball class
Initialises Simulation with passed initial conditions
Runs for 2000 frames to thermalise
Reset total_time and total_mom so no data is taken up to this point
Run for 500 more frames and then record pressure

"""

from thermoCheck import *
import numpy as np 
import matplotlib.pyplot as plt 
from tqdm import trange

def pt(cont,n,m,r,v_bottom,v_top):
    '''
    See module docstring for information on purpose
    
    Parameters:
        cont = Container radius
        n = number of balls
        m - mass of balls
        r = radius of balls
        v_bottom = lower bound of velocity component speed
        v_top = upper bound of velocity componenet speed
    
    
    Momentum data is erased after 2000 frames 
    Simulation run for 500 more frames
    Pressure is only dependent on final 500 frames
    
    Returns temperature and pressure arrays to be plotted against each other
    '''
    assert isinstance(n,int)
    assert isinstance(cont,float)    
    assert isinstance(m,float)
    assert isinstance(r,float)
    s=Sim(cont)
    temp_array=[]
    pres_array=[]
    
    for i in trange(1,21):
        s.generateBalls(n,m,r, 10*i + v_bottom, 10*i + v_top) 
        s.Run(2000)

        s.total_time=0. # reset time and momentum 
        s.total_mom=0.
        s.Run(500) # collect pressur edata over next 500 collisions
        t=Temp(s)
        p=Pres(s)
        temp_array.append(t)
        pres_array.append(p)

    return temp_array, pres_array
    

    
    