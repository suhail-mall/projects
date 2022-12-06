# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 00:33:21 2019

@author: suhai

Contains only Sim class
Run help() on Sim class for docstring

"""
from initBall import *
import numpy as np
import pylab as pl
from tqdm import trange
 

class Sim:
    '''
    Simulation class is initialised with a radius which is used to initialise a Container (derived calss of Ball).
    (n) Balls can then be generated with associated mass and radius. Velocity is random chosen using a normal distribution and position is systematic
    When the Run() method is called: next collision is found, time advanced to that point, collision carried out.
    Run() takes frame number (N) --> N collisions are carried out.
    
    Note that tqdm package must be installed 
    -------------------------------------------------------------------------------------------------------------------------------------------------------

    Example of use in iPython terminal:
        
        s=Sim(10.)
        s.generateBalls(100,1.,0.1,-10.,10.)
        s.Run(2000,True)
    --------------------------------------------------------------------------------------------------------------------------------------------------------     
        
    '''    
    def __init__(self, rad):
        '''
        Initialises Simulation object with a Container object.
        Note that passed Container radius value 'rad' must be type float.
        Initialises all variables as floats
        
        mass and radius of balls is property of Simulation (all the same) so stored here
        total_time and total_mom are used to calculate pressure exerted on Container. Property of simulation so stored here.
        '''
        assert isinstance(rad,float)
        self.container_rad          = rad
        self.container              = Ball(1E12,rad,[0.,0.],[0.,0.])
        self.container.is_container = True
        self.balls      = [] # initalise, will hold all balls in Simulation object
        self.num_balls  = 0
        self.ball_mass  = 0.
        self.total_time = 0.
        self.total_mom  = 0.
    
    def generateBalls(self, n, m, r, v_bottom, v_top):
        '''
        Used to initialise a large number of Ball objects within a Simulation.
        Note that all passed variables must be type float.
        -------------------------------------------------------------------------------------------------------------------------------------------------------
        
        n= number of balls
        
        m= mass of each ball
        
        r= radius of each ball
        
        v_bottom= lower bound for velocity component speed
        
        v_top= upper bound for velocity component speed
        -------------------------------------------------------------------------------------------------------------------------------------------------------
        
        Can:
            Set spacing of Balls.
            Change distribution from which ball velocity components are chosen, e.g. uniform, normal, etc
        '''
        spacing  = 10*r # spacing between balls. Should reduce for large ball_radius or large ball number
        cont_rad = self.container_rad
        self.num_balls  = n
        self.ball_mass  = m
        self.total_mom  = 0. # reset these each time a new set of balls is generated
        self.total_time = 0.
        
        ball_array = [self.container] # to hold all Ball objects
                
        grid_list  = np.arange(-cont_rad+0.1, cont_rad-0.1, spacing)
        final_grid = []
        
        for i in grid_list:
            for j in grid_list:
                mag = np.sqrt( (i)**2 + (j)**2 ) +r +0.01 
                if mag<cont_rad:
                    final_grid.append([ float(i),float(j) ])
        if len(final_grid)<n:
            raise Exception('too many balls for container size')
            
        for k in range(0,n):
            gen_vx  = np.random.normal(v_bottom,v_top)
            gen_vy  = np.random.normal(v_bottom,v_top)
            gen_vel = [gen_vx, gen_vy]
            
            gen_pos = final_grid[k]
            
            ball_array.append(Ball(m,r,gen_pos,gen_vel))
        
        self.balls = ball_array
        return ball_array
    
    
    def lowestTime(self): 
        '''
        Iterates through each non-repeating combination of two balls and uses the initBall.timeToColl() method to calculate the time to collsion.
        Records the time to the first collision as well as the two balls colliding in that collision.
        Lowest time and balls are output 
        Returns lowest_time and balls involved in next collision
        -------------------------------------------------------------------------------------------------------------------------------------------------------  
        '''
        N=len(self.balls)
        #use coll_1 to denote ball_1 of first collision
        #use coll_2 to denote ball_2 of first collision
        
        lowest_time=1E20 # this is dodgy but need an initial value to check against
        for i in range(N): 
            for j in range(i+1,N): # to prevent unnecessary calculations - see lab book
                t, is_coll = self.balls[i].timeToColl(self.balls[j])
               
                if is_coll and t <= lowest_time: # Collision has happened and is first collision
                    lowest_time = t
                    coll_1 = i
                    coll_2 = j    
                   
        return lowest_time, coll_1,coll_2         
 
    def nextCollision(self):
        '''
        Moves all balls to next collsion.
        Collision carried out.
        Simulation advanced by small time to prevent sticking
        Does not return any object
        
        If ball_1 is object Container, momentum change of ball_2 is recorded to calculate pressure exerted on Container.
        This calculation is done in module thermoCalc
        -------------------------------------------------------------------------------------------------------------------------------------------------------       
        '''        
        N         = len(self.balls)
        t_to_coll = self.lowestTime()[0]
        coll_1    = self.lowestTime()[1]
        coll_2    = self.lowestTime()[2]
        
        self.total_time += t_to_coll
        
        if coll_1 == 0:
            v_before = self.balls[coll_2]._vel.copy() # save the velocity of the ball before collision so that we can calculate the change in momentum
        
        [ball.move(t_to_coll) for ball in self.balls]
        
        self.balls[coll_1].collide(self.balls[coll_2])
        
        if coll_1==0: # If involve container, use to calculate pressure 
            v_after = self.balls[coll_2]._vel.copy()
            m       = self.ball_mass
            dv      = v_after - v_before
            dv_mag  = np.sqrt( (dv[0])**2 + (dv[1])**2 )
            dv_mom  = m*dv_mag
            self.total_mom += dv_mom
            
        [ball.move(1e-12) for ball in self.balls]
           
    def drawFrame(self):
        '''
        Draws a single frame.
        Adds a container and then iterates through balls, adding circular patch of their radius at their position.
        Container patch is not filled
        First ball in ball_array is coloured blue to help debugging 
        -------------------------------------------------------------------------------------------------------------------------------------------------------
        '''
        N = len(self.balls)
        f = pl.figure()
        
        cont_rad = self.container_rad
        ax = pl.axes(xlim = (-1*cont_rad,cont_rad), ylim =(-1*cont_rad,cont_rad)) 
        pl.gca().set_aspect('equal',adjustable='box')
        
        patch_container = pl.Circle([0.,0.], cont_rad, ec='b', fill=False, ls='solid')
        ax.add_patch(patch_container)
        
        patch = pl.Circle(self.balls[1]._pos, self.balls[1].rad, fc='b')
        ax.add_patch(patch)
        
        for i in range(2,N):
            patch = pl.Circle(self.balls[i]._pos, self.balls[i].rad, fc='r')
            ax.add_patch(patch)           
        f.show() 

        
    def Run(self,no_frames,animate = False):
        '''
        Draw frame and carry out next collision.
        If animate ==True, will continue this for specified number of frames
        Frames update to next collision.
        -------------------------------------------------------------------------------------------------------------------------------------------------------
        
        no_frames= number of collisions to run simulation for. Must be type integer
        animate= will animate the simulation if True. Must be type Boolean
        
        '''

        if animate: #if True
            self.drawFrame()
        for frame in trange(no_frames): # trange can be used here to check progress of simulation 
            self.nextCollision()
            if animate:
                pl.pause(0.0001)
        if animate:
            pl.show()

 