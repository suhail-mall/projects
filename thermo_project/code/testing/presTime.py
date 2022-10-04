# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 19:48:52 2019
@author: ssm2617

Takes initialisation parameters for a simulation and initlises simulation.
Samples the pressure exerted on the system every 500 total collisions (not just collisions on the container wall).
THIS IS CUMULATIVE, I.E. THE TOTAL TIME AND MOMENTUM IMPARTED ON THE WALL DOES NOT RESET FOR EACH SAMPLE. 
"""

from thermoCalc import *
import numpy as np 
import matplotlib.pyplot as plt 
from tqdm import trange

def presFrame(cont,n,m,r,v_bottom,v_top):
    '''
    Takes initialisation parameters for a simulation and initlises simulation.
    Samples the pressure exerted on the system every 500 total collisions (not just collisions on the container wall).
    THIS IS CUMULATIVE, I.E. THE TOTAL TIME AND MOMENTUM IMPARTED ON THE WALL DOES NOT RESET FOR EACH SAMPLE. 
    
    Takes parameters:
        cont = container radius
        n = number of balls
        m = mass of balls
        r = radius of balls
        v_bottom = lower bound of ball velocity-component
        v_top = upper bound of ball velocity-component        
        
    Returns lists that can be plotted against each other
    '''
    assert isinstance(n,int)
    assert isinstance(cont,float)    
    assert isinstance(m,float)
    assert isinstance(r,float)
    s=Sim(cont)
    pres_array=[]
    frame_array=[]
    
    s.generateBalls(n,m,r,v_bottom,v_top)
    
    for i in trange(1,21):        
        s.Run(500)
        t=i*500
        p=Pres(s)
        frame_array.append(t)
        pres_array.append(p)
    
    return frame_array, pres_array