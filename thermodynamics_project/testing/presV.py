# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 12:49:58 2019

@author: ssm2617

Passed parameters required to initialise Sim and Ball class
Initialises Simulation with passed initial conditions
Runs for 2000 frames to thermalise
Measure pressure

Increase container radius by value 1
Run for 20 data points
"""

from thermoCalc import *
import numpy as np 
import matplotlib.pyplot as plt 
from tqdm import trange

def pV(cont,n,m,r,v_bottom,v_top):
    assert isinstance(n,int)
    assert isinstance(cont,float)    
    assert isinstance(m,float)
    assert isinstance(r,float)
    
    pres_array=[]
    v_array=[]
    
    for i in trange(0,21):
        cont_var=cont+i
        s=Sim(cont_var)
        s.generateBalls(n,m,r,v_bottom,v_top)
        s.Run(3000)
        v_var = np.pi * (cont_var)**2
        v_array.append(v_var)
        p=Pres(s)
        pres_array.append(p)
    return v_array, pres_array