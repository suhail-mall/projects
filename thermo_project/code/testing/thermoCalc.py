# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 19:31:34 2019
@author: ssm2617

Contains methods that take a Sim class and calculate the kinetic energy of Balls in system
Can relate this to temperature of system

"""
import numpy as np 
from nBall import *
from tqdm import trange

def Ek(sim):
    '''
    Is passed a Sim class
    Iterates through each ball in the Sim class noting its mass and velocity. 
    Uses these to calculate total kinetic energy
    '''
    ek=0
    N=sim.num_balls
    for i in range(0,N):
        vel=sim.balls[i]._vel
        speed_sq=(vel[0])**2 + (vel[1])**2
        mass=sim.balls[i].mass
        ek+=0.5*mass*speed_sq
    return ek

def Temp(sim):
    '''
    Calls Ek function to calculate energy of simulation
    Relates this to temperature of system
    '''
    N=sim.num_balls
    t=(2)*Ek(sim)/(N*1.38064852E-23) 
    return t
   
def Pres(sim):
    '''
    Passed Sim class
    Accesses total_time and total_mom from Sim
    Returns pressure
    '''
    return sim.total_mom / (sim.total_time* 2*np.pi *sim.container_rad)

     
def ekTest(sim):
    '''
    Quick test to ensure Kinetic Energy is being conserved
    Ek should not change between collisions
    
    Passed Sim class
    Advances simulation by 10 frames
    calculate Ek of system
    Stores Ek for frame 
    Repeat for 500 data points
    '''
    ek_array=[]
    for i in trange(0,500):
        sim.Run(10)
        ek=Ek(sim)
        ek_array.append(ek)
    return ek_array