# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:32:27 2019

@author: suhailMall
"""
import numpy as np
import matplotlib.pyplot as plt 
import urllib as url 
import math 
from tqdm import *
import scipy.optimize as op


#%%    READ DATA FROM ONLINE TEXT FILE (REQUIRES IMPERIAL LOGIN)
# username = 'ssm2617'
# file = url.request.urlopen('https://www.hep.ph.ic.ac.uk/~ms2609/CompPhys/neutrino_data/%s.txt'%username) # Request data 

# data = [x for x in file] # function returns object, not list
# N = int ((len(data) - 5)/2) # spaces and labels in text file 

# k_arr = [float(str(k)[2:-3]) for k in data[2:N+2]] # Based on webpage format
# lambda_arr = [float(str(r)[2:-3]) for r in data[N+5:]]
    
#%%        READ DATA FROM TEXT FILE (PRE-DOWNLOADED AND WRITTEN OUT)    
k_arr = np.loadtxt('data/k_data.txt')
lambda_arr = np.loadtxt('data/rate_data.txt')

#%%        OSCILLATION PROBABILITY AND LIKELIHOOD 
#          FOR 2D PARAM SPACE (theta_{23} AND m_squared_{23}) 
           
def OscProb(E, theta, m, L=295):
    ''''
    > Probability of a neutrino oscillating from a \mu neutrino to a 
    \mu neutrino
    > For neutrino of energy E and oscialltion parameters theta and m_squared 
    -------
    PARAMS 
    E : energy of neutrino 
    theta : mixing angle of neutrinos 2 and 3 
    m : mass-squared difference of neutrinos 2 and 3 
    L : length of detector, taken to be 295m 
    -------
    RETURNS 
    P_11 : probability of neutrino oscillating from flavour 1 to flavour 1 
    '''
    return 1 - (np.sin(2*theta))**2 * (np.sin(1.267*m*L/E))**2    
 
def poisson(mean, obs):
    '''
    > Poisson Distribution 
    '''
    return mean**obs * np.exp(-mean)/math.factorial(int(obs))
    
    
def lambda_i_2d(lambda_pred, E, theta, m, L=295): 
    '''
    > Multiply the predicted observations (with no oscillation) by the probability of 
    oscillation to the desired neutrino flavour 
    > lmabda_i because it will be used as the observations for the i'th energy bin 
    > this model depends on two parameters: theta and m 
    -------
    PARAMS 
    lambda_pred : predictedobservations without neutrino oscillation 
    E : energy of neutrino 
    theta : mixing angle of neutrinos 2 and 3 
    m : mass-squared difference of neutrinos 2 and 3 
    L : length of detector, taken to be 295m 
    -------
    RETURNS 
    lambda_i : lambda of energy bin i, with neutrino oscillations 
    '''
    return lambda_pred * OscProb(E, theta,m)  


def NLL_2d(theta,m):
    '''
    > Negative Log Likelihood of neutrino flavour oscillation modulating 
    the predicted observations as described by OscProb()
    > This model only depends on the mixing angle and mass-squared difference
    of neutrinos 2 and 3
    -------
    PARAMS 
    theta : mixing angle of neutrinos 2 and 3 
    m : mass-squared difference of neutrinos 2 and 3 
    -------
    RETURNS 
    nll : the negative log likelihood of these parameters describing the data 
    '''
    E_arr = np.linspace(0.025,9.975,200)# bin energy is centre of range 
    nll = -1*np.sum(
        [np.log(poisson(lambda_i_2d(lambda_arr[i], E_arr[i], theta, m), k_arr[i]))
                   for i in range(200)])
    return nll


#%%        NLL FOR 3D PARAM SPACE (theta_{23}, m_squared_{23}, AND cs)
#          cs multiplies the oscillation probability with a term linear in E

def lambda_i_3d(lambda_pred, E, theta, m, cs, L=295): 
    '''
    > Multiply the predicted observations (with no oscillation) by the probability of 
    oscillation to the desired neutrino flavour 
    > lmabda_i because it will be used as the observations for the i'th energy bin 
    > this model depends on 3 parameters: theta, m, and cs 
    -------
    PARAMS 
    lambda_pred : predicted observations without neutrino oscillation 
    E : energy of neutrino 
    theta : mixing angle of neutrinos 2 and 3 
    m : mass-squared difference of neutrinos 2 and 3 
    cs : a constant of proportionality for the oscilaltion probability and 
    energy 
    L : length of detector, taken to be 295m 
    -------
    RETURNS 
    lambda_i : lambda of energy bin i, with neutrino oscillations 
    '''
    return lambda_i_2d(lambda_pred, E, theta, m) *  (cs * E)


def NLL_3d(theta,m,cs):
    '''
    > Negative Log Likelihood of neutrino flavour oscillation modulating 
    the predicted observations as described by OscProb()
    > This model depends on the mixing angle and mass-squared difference
    of neutrinos 2 and 3, as well as a term proportional with the energy with
    a constant factor cs
    -------
    PARAMS 
    theta : mixing angle of neutrinos 2 and 3 
    m : mass-squared difference of neutrinos 2 and 3 
    cs : a constant of proportionality for the oscilaltion probability and 
    energy 
    -------
    RETURNS 
    nll : the negative log likelihood of these 3 parameters describing the data 
    '''
    E_arr = np.linspace(0.025,9.975,200)# bin energy is centre of range 
    nll = -1*np.sum(
        [np.log(poisson(lambda_i_3d(lambda_arr[i], E_arr[i], theta, m, cs), k_arr[i]))
                   for i in range(200)])
    return nll

#%%        HESSIAN and COV MATRICES 

def Hessian(f,x,eps):
    '''
    > Calculates the second derivative of an n-dim function at position x
    using a central difference scheme, given as the Hessian matrix 
    H_{ij} = \partial_i \partial_j f(x)
    -------
    PARAMS 
    f : a function that takes any number of arguments 
    x : the position in parameter space where the derivative is calculated 
    eps : an n-dim vector for finite differences along each parameter 
    -------
    RETURNS 
    H : Hessian matrix of size (n,n)
    '''
    eps_mat = np.diag(eps)
    n = len(x)
    H = np.zeros([n,n])
    for i in range(n):
        for j in range(n): # Calculated on paper 
            H[i,j] = 1/(4*eps[i]*eps[j]) \
               * (f(*(x+eps_mat[i]+eps_mat[j])) \
               -  f(*(x+eps_mat[i]-eps_mat[j])) \
               -  f(*(x-eps_mat[i]+eps_mat[j])) \
               +  f(*(x-eps_mat[i]-eps_mat[j]))) 
    return H

def cov(H): 
    '''
    > Covariance matrix taken as the inverse of half the Hessian matrix 
    -------
    PARAMS 
    H : the Hessian matrix 
    -------
    RETURNS 
    c : the covariance matrix 
    '''
    return np.linalg.inv(0.5*H)

#%%        Plot Oscillation Probability and Expected vs Observed Events 
def PlotOscProbE(theta = np.pi/4,m_squared=2.4e-3,L=295, plot = True):
    E = np.linspace(0,10,200)
    prob = [OscProb(x,L,theta,m_squared) for x in E]
    if plot:
        plt.plot(E, prob, '-')
        plt.grid()
        plt.xlabel('E, GeV')
        plt.ylabel('$P_{osc} (E)$')
   
def Hists():
    plt.clf()
    E = np.linspace(0,10,200)
    
    plt.figure(1)
    plt.grid()
    plt.xlabel('E, GeV')
    plt.ylabel('Number Observed Events and observations With Oscillation')
    plt.bar(E,k_data,width=0.05, label = 'Number of Observed Events')

    OR = [lambda_i(lambda_arr[i],E[i], 295, np.pi/4, 2.4e-3) for i in range(200)]
    plt.bar(E,OR, width = 0.05, color = 'orange', label = 'Simulated observations with oscillation')
    plt.legend()
    
    plt.figure(2)
    plt.grid()
    plt.xlabel('E, GeV')
    plt.ylabel(' observations Without Oscillation (Simulation)')
    plt.bar(E, lambda_arr, width = 0.05)
    
def HeatMap(length = 500, plot = True):
    theta_arr = np.linspace(0,np.pi/2,length)
    m_arr = np.linspace(0,5e-3,length)
    nll_arr = np.zeros([length,length])
    
    for i in trange(length):
        for j in range(length):
            nll_arr[i,j] = NLL_2d(theta_arr[i],m_arr[j])
        
    if plot:
        fig, ax = plt.subplots()
        plt.imshow(nll_arr, extent = [0,5e-3,0,np.pi/2], aspect = 'auto')
    return fig, ax 


