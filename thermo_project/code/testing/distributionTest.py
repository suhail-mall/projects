# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:07:12 2019

@author: ssm2617

Calculates:
    Distance of each ball from centre of Simulation container
    Speed of each ball in Simulation

Returns array that can be used to plot Histogram of each distribution

Run simulation for largne number of collisions then test distributions

Should produce: distance distribution ~ r
    Speed is Maxwell-Boltzmann
"""
import numpy as np
import matplotlib.pyplot as plt


from nBall import *
from tqdm import trange

def testDistribution(Sim):
    '''
    Passed a Sim class, looks at the final position and velocity arrays of each ball and constructs an array of each. A histogram can then be plotted
    '''
    dist_list=[]
    speed_list=[]
    N=Sim.num_balls
    for i in range(N):
        ball_pos=Sim.balls[i]._pos        
        dist=np.sqrt((ball_pos[0])**2 + (ball_pos[1])**2)
        dist_list.append(dist)
        
        ball_vel=Sim.balls[i]._vel
        speed=np.sqrt((ball_vel[0])**2 + (ball_vel[1])**2)
        speed_list.append(speed)
    return dist_list, speed_list  

def Hist(array,Bins=50):
    '''
    Shows a histogram after being passed an array
    For convenience of not having to type plt
    '''
    plt.plot()
    plt.hist(array, bins=Bins) 