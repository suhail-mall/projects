# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 21:21:59 2019

@author: ssm2617


Used to initialise Ball class and Container derived class. 
Contains methods to check time to next collision and update velocity due to collision.
Use help() on Ball class for docstring
""" 
import numpy as np


class Ball:
    '''
    Initialises an object with 
	mass
	radius
	position
	velocity.
    Contains methods for Ball object
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Can be used directly by user but mostly accessed by Sim class in nBall module
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    '''
    def __init__(self, mass, rad, pos, vel):
        '''
        Initialises Ball object with variables:
            mass
            radius
            position array (shape 1x2, element type float)
            velocity array (shape 1x2, element type float)
            is_container = False (can't initalise Ball as Container derived class)
        '''
        import numpy as np
        self.mass = float(mass)
        self.rad = float(rad)
        self.is_container = False # Have derived class Container - sets to True
        
        self.check(pos,vel)
        self._pos = np.array(pos, dtype = float)
        self._vel = np.array(vel, dtype = float)       
    
    def check(self,pos_array, vel_array):
        '''
        Used to check the length and type of element in pos and vel arrays.
        
        '''
        if len(pos_array) !=2:
            raise Exception('position array must have length 2')    #array length check
        if len(vel_array) !=2:
            raise Exception('velocity array must have length 2')
        for i in range(0,2):
            if type(pos_array[i]) != float:
              raise Exception('position array elements must be of type float')
            if type(vel_array[i]) != float:
                raise Exception('velocity array elements must be of type float')
    
    def showPos(self):
        '''
        Returns position array of object
        '''
        return self._pos
    
    def showVel(self):
        '''
        Returns velocity array of object
        '''
        return self._vel
    
    def move(self,dt):
        '''
        Moves the ball to a time in the future determined by user
        Simple calculation because Ball object does not experience acceleration if not colliding
        '''
        # self._pos = np.array([(pos + vel*dt) for pos,vel in zip(self._pos, self._vel)]) # Wanted to do with zip funcn 
        self._pos = self._pos + self._vel*dt
    
 
    def timeToColl (self, other):
        '''
        Calculates the time to collide based on the positions and velocities of two objects
        Takes other ball as input 
        
        If time to collision is valid (i.e. real and positive), return time and True boolean
        If time to collision is invalid, return 0. and False boolean
        Will only consider time if passed Boolean is True
        '''
        import numpy as np
        r = self._pos - other._pos
        v = self._vel - other._vel
        #if one of the balls is a container, then the formula for time to collide is different 
        if self.is_container:
            R = self.rad - other.rad
        elif other.is_container:
            R = other.rad - self.rad
        else:
            R = self.rad + other.rad   
        r_dot_v  = np.dot(r,v)
        r_mag_sq = np.dot(r,r)
        v_mag_sq = np.dot(v,v)
        
        is_coll = True # assume going to collide until shown not to
        dt      = 0.
        #dt gives two solns --> see lab book for details of cases
        root = (r_dot_v)**2 -(v_mag_sq)*(r_mag_sq - R**2)
        if root>0: #collision time is real 
            dt_pos =  (- r_dot_v + np.sqrt(root)) / (v_mag_sq)
            dt_neg =  (- r_dot_v - np.sqrt(root) )/ (v_mag_sq)
            #If both positive, take smallest positive 
            #If one pos and one neg, have one is container --> choose pos soln
            if (dt_pos>0 and dt_neg>0): #both dt are positive --> take smallest
                if dt_pos < dt_neg:
                    dt = dt_pos
                elif dt_neg < dt_pos:
                        dt = dt_neg
                else:
                    dt=dt_pos # if both are equal
            elif (dt_pos>0 and dt_neg<0) or (dt_pos<0 and dt_neg>0): # one is a container --> take positive soln
                if dt_pos > 0: 
                    dt=dt_pos
                elif dt_neg>0:
                    dt=dt_neg
            elif dt_pos == dt_neg: # just in case - SHOULDN'T BE HERE
                dt = dt_pos 
            elif dt_pos<0 and dt_neg<0:#both solns negative
                is_coll =False
        else: # if soln isn't real (root<0) 
            dt = 0.
            is_coll = False
        return dt, is_coll # still need to pass dt as a float
    
    
    def normalise(r):
        '''
        Normalise a passed vector (will only be passed 1x2 arrays by code)
        '''
        r = np.array(r) # Ensure use of numpy array operations 
        return r/np.sqrt(np.dot(r,r)) 
        
    def collide(self,other):  
        '''
        Updates the velocity of two objects after a collision has occured.
        Takes other Ball object
        Return both objects
        '''       
        u1 = self._vel
        r1 = self._pos
        u2 = other._vel
        r2 = other._pos
        
        sum_mass = self.mass + other.mass
        u12 = u1 - u2
        r12 = r1 - r2
        
        self._vel  -= (2*other.mass/sum_mass) * np.dot(u12,r12) * (r12)  / (np.dot(r12,r12))       
        other._vel -= (2*self.mass/sum_mass)  * np.dot(u12,r12) * (-r12) / (np.dot(r12,r12))  
        return self, other


