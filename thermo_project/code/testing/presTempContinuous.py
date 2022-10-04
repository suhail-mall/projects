# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 20:20:58 2019

@author: ssm2617

Passed parameters required to initialise Sim and Ball class
Initialises Simulation with passed initial conditions
Runs for 2000 frames to thermalise
Records cumulative pressure over this time

Initalises simulations with initial velocity componenet speed distribution over different range
Here translates range by +10 for each simulation
"""

from thermoCalc import *
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
    
    
    Pressure measured here is cumulative, i.e. measured continuously from start of simulation 
    
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
        s.generateBalls(n,m,r, 10*i + v_bottom, 10*i + v_top) #generate new balls. This also resets total time and total momentum        
        s.Run(2000)
        t=Temp(s)
        p=Pres(s)
        temp_array.append(t)
        pres_array.append(p)
 
    return temp_array, pres_array       
