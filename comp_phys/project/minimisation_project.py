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
username = 'ssm2617'
file = url.request.urlopen('https://www.hep.ph.ic.ac.uk/~ms2609/CompPhys/neutrino_data/%s.txt'%username) # Request data 

data = [x for x in file] # function returns object, not list
N = int ((len(data) - 5)/2) # spaces and labels in text file 

k_arr = [float(str(k)[2:-3]) for k in data[2:N+2]] # Based on webpage format
lambda_arr = [float(str(r)[2:-3]) for r in data[N+5:]]
    
#%%        READ DATA FROM TEXT FILE (PRE-DOWNLOADED AND WRITTEN OUT)    
k_arr = np.loadtxt('data/k_data.txt')
lambda_arr = np.loadtxt('data/rate_data.txt')

#%%        OSCILLATION PROBABILITY AND LIKELIHOOD 
def OscProb(E, L,theta,m_squared):
    ''''
    > Probability of a neutrino oscillating from a \mu neutrino to a 
    \mu neutrino
    > For neutrino of energy E and oscialltion parameters theta and m_squared 
    
    '''
    return 1 - (np.sin(2*theta))**2 * (np.sin(1.267*m_squared*L/E))**2    
   
def lambda_i(lambda_pred, E, theta, m,L=295): 
    '''
    > Multiply the predicted rate (with no oscillation) by the probability of 
    oscillation to the desired neutrino flavour 
    > lmabda_i because it will be used as the rate for the i'th energy bin 
    '''
    return lambda_pred * OscProb(E, L, theta,m)  

def poisson(lambda_pred, E, theta, m, k):
    mean = lambda_i(lambda_pred,E,theta,m)  
    return mean**k * np.exp(-mean)/math.factorial(int(k))

def NLL(theta,m):
    E_arr = np.linspace(0.025,9.975,200)# bin energy is centre of range 
    nll = -1*np.sum([np.log(poisson(lambda_arr[i],E_arr[i],theta,m,k_arr[i]))
                   for i in range(200)])
    return nll
    
    

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
    plt.ylabel('Number Observed Events and Rate With Oscillation')
    plt.bar(E,k_data,width=0.05, label = 'Number of Observed Events')

    OR = [lambda_i(rate_pred[i],E[i], 295, np.pi/4, 2.4e-3) for i in range(200)]
    plt.bar(E,OR, width = 0.05, color = 'orange', label = 'Simulated Rate with oscillation')
    plt.legend()
    
    plt.figure(2)
    plt.grid()
    plt.xlabel('E, GeV')
    plt.ylabel(' Rate Without Oscillation (Simulation)')
    plt.bar(E, rate_pred, width = 0.05)
    
def HeatMap(length = 500, plot = False):
    theta_arr = np.linspace(0,np.pi/2,length)
    m_arr = np.linspace(0,5e-3,length)
    nll_arr = np.zeros([length,length])
    
    for i in trange(length):
        for j in range(length):
            nll_arr[i,j] = NLL(theta_arr[i],m_arr[j])
        
    if plot:
        plt.figure()
        plt.imshow(nll_arr, extent = [0,5e-3,0,np.pi/2], aspect = 'auto')
    return theta_arr, m_arr, nll_arr 
    
    
#%%    1D MINIMISATION 
    
def MinPoint(x0,x1,x2, y0,y1,y2): 
    '''
    > Given three data points, fit to a quadratic and find the minimum
    '''        
    return 0.5 * ( (x2**2 - x1**2)*y0 + (x0**2 - x2**2)*y1 + (x1**2 - x0**2)*y2) / ((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2)


def ParaMinTheta(guess_points = [0.8,0.9,0.95], m=2.4e-3):
    '''
    > Parabolic minimisation of NLL varying just the theta parameter 
    > Approximate as parabola and find minimum point, use lowest three points 
    to fit the next parabola. Repeat until finishing condition met 
    '''
    change = 1. # initialised so that it doesn't terminate immediately 
    x0,x1,x2 = np.sort(guess_points)
    xmin_old = x0
    
    while change > 1e-3:           
        xmin =  MinPoint(x0,x1,x2,NLL(x0,m),NLL(x1,m),NLL(x2,m))

        points = [x0,x1,x2,xmin]
        nll_points = [NLL(p,m) for p in points]
        
        for point in points:
            if np.max(nll_points) == NLL(point,m):
                points.remove(point)
                break 
        x0,x1,x2 = np.sort(points) 
        change = np.abs(xmin_old - xmin)/xmin_old
        xmin_old = xmin
    return xmin,[x0,x1,x2]  # x_min doesn;t necessarily belong to [x0,x1,x2]

def ParaMinM(guess_points = [2e-2,2.5e-3,3e-3], theta = 0.95):
    '''
    > Parabolic minimisation of NLL varying just the m parameter 
    > Approximate as parabola and find minimum point, use lowest three points 
    to fit the next parabola. Repeat until finishing condition met 
    '''
    change = 1.
    x0,x1,x2 = np.sort(guess_points)
    xmin_old = x0
    
    while change > 1e-3:           
        xmin =  MinPoint(x0,x1,x2,NLL(theta,x0),NLL(theta,x1),NLL(theta,x2))

        points = [x0,x1,x2,xmin]
        nll_points = [NLL(theta,p) for p in points]
        
        for point in points:
            if np.max(nll_points) == NLL(theta,point):
                points.remove(point)
                break    
        x0,x1,x2 = np.sort(points) 
        change = np.abs(xmin_old - xmin)/xmin_old
        xmin_old = xmin
    return xmin,[x0,x1,x2]  


#%%        Uncertainty 

def UncGauss(theta_min, sig_guess=0.05):
    '''
    > Treat likelihood as Gaussian, then \sigma is approximately found as the 
    value where the log-likelihood decreases by an absolute value of 0.5
    > Seek and return the values either side of the minimum point that satisfy
    this using scipy.optimize.fsolve
    '''
    nll_min = NLL(theta_min, 2.4e-3)
    theta_minus = op.fsolve(lambda theta : NLL(theta,2.4e-3)-nll_min-0.5 , 
              theta_min-sig_guess, full_output=True)
    theta_plus  = op.fsolve(lambda theta : NLL(theta,2.4e-3)-nll_min-0.5 , 
              theta_min+sig_guess, full_output=True)
    
    sig_minus = theta_min - theta_minus
    sig_plus  = theta_plus - theta_min
    return sig_minus, sig_plus 


def UncDeriv(t_min,h = 1e-7):
    '''
    > 
    '''
    m = 2.4e-3
    deriv2_t = (NLL(t_min + h, m) - 2*NLL(t_min, m) + NLL(t_min - h, m) ) / h**2
    return np.sqrt(1 / deriv2_t)

def UncPara(t_min,m,t0,t1,t2):
    y0 = NLL(t0,m)
    y1 = NLL(t1,m)
    y2 = NLL(t2,m)    
    return (2*(y0/((t0-t1)*(t0-t2)) + y1/((t1-t0)*(t1-t2)) + y2/((t2-t0)*(t2-t1))))**-0.5

#%%    UNIVARIATE 

def Univariate(guess_theta = [0.8,0.9,0.95],guess_m=[2e-3,2.5e-3,3e-3]):
    finished = False
    t_points = np.sort(guess_theta) 
    m_points = np.sort(guess_m)
    nll_old = 1e11
    t_min = t_points[1]
    m_min = m_points[1]
    t_tests = [t_min]
    m_tests = [m_min]
    
    while not finished:

        m_min, m_points = ParaMinM(m_points, t_min)
        t_tests.append(t_min)
        m_tests.append(m_min)        
        t_min, t_points = ParaMinTheta(t_points, m_min)
        t_tests.append(t_min)
        m_tests.append(m_min)            
        
        nll = NLL(t_min,m_min)
        ratio = 1 - nll/nll_old
        nll_old = nll      
        if  np.abs(ratio) < 5e-8:
            finished = True 
    return [t_min, m_min], [t_tests, m_tests]

       

def UncDeriv2(t_min,m_min,h=1e-7):
    deriv2_t = (NLL(t_min + h, m_min) - 2*NLL(t_min, m_min) + NLL(t_min - h, m_min) ) / h**2
    deriv2_m = (NLL(t_min, m_min + h) - 2*NLL(t_min, m_min) + NLL(t_min, m_min - h) ) / h**2
    return np.sqrt(1 / deriv2_t), np.sqrt(1 / deriv2_m)
#%%        Gradient Minimisation 
    
def GradMin(test_t = 0.8, test_m = 3e-3,alphas = [1e-4,1e-9]):
    finished = False 
    eps = 1e-7
    alpha_t, alpha_m = alphas
    old_t = 1e11 # initialise
    old_m = 1e11
    t_arr = [test_t]
    m_arr = [test_m]
    
    while not finished:        
        grad_t = 0.5 * (NLL(test_t + eps, test_m) - NLL(test_t - eps, test_m)) / eps # calculate gradient
        grad_m = 0.5 * (NLL(test_t, test_m + eps*1e-3) - NLL(test_t, test_m - eps*1e-3)) / (eps*1e-3)
               
        test_t -= alpha_t*grad_t # step away 
        test_m -= alpha_m*grad_m
        
        t_arr.append(test_t)
        m_arr.append(test_m)
        
        ratio = 1 - NLL(test_t, test_m)/NLL(old_t,old_m)       
        if 0 < ratio < 1e-8:
            finished = True 
        else:
            old_t = test_t
            old_m = test_m
    return [test_t, test_m], [t_arr,m_arr]

#%%    3d minimisation 

def LambdaI3(lambda_pred, E, L, theta, m, cs ): # the predicted lambda * osc prob
    return OscProb(E, L, theta,m) * lambda_pred * (cs * E)

def NLL3(theta,m,cs): # Updated to take cross-section term 
    E_array = np.linspace(0.025,9.975,200)# bin energy is centre of range
    nll = 0. # initialise NLL element 
    for i in range(200): # sum up terms 
        lambda_i = LambdaI3(rate_pred[i],E_array[i], 295, theta, m,cs) # fixed L
        nll += -np.log((lambda_i**k_data[i]) * np.exp(-1*lambda_i) / math.factorial(k_data[i]) )
    return nll
    
def QN(f = NLL3 ,x0 = [0.86,2.88e-3,1.], alpha = [1e-4,1e-9, 1e-3]):
#    alpha = [1e-9,1e-9]
#    alpha = [1e-9,1e-9]
    n = len(x0) # number of dimensions 
    x0 = np.array(x0)
    g0 = np.identity(n)
    finished = False
        
    x_arr = []  
    grad_arr = []
    gamma_arr = []
    delta_arr = []
    g_arr = []
    
    if n == 2:
        grad0 = Grad2(f,x0)
    if n == 3:
        grad0 = Grad3(f,x0)
    
    grad_arr.append(grad0) # easy to access in loops instead of tracking current and old 
    x_arr.append(x0)
    g_arr.append(g0)   
    
    while not finished: #while it's not finished
        # Update position
        
        a =  alpha *grad_arr[-1]
        x_arr.append( x_arr[-1] -  np.matmul(g_arr[-1],a) )  # calculate x_{n+1} then calculate quantities at x_{n+1}

        # Update other quantities 
        if n==2:
            grad_arr.append(Grad2(f,x_arr[-1])) # evaluated at x_{n+1}
        if n==3:
            grad_arr.append(Grad3(f,x_arr[-1]))
            
        gamma_arr.append(grad_arr[-1] - grad_arr[-2])
        delta_arr.append(x_arr[-1] - x_arr[-2])
        
        g_arr.append(Gnext(g_arr[-1], delta_arr[-1], gamma_arr[-1]))
                
        if n ==2:  
            x0_new = x_arr[-1][0] # x_i for eachdirection 
            x1_new = x_arr[-1][1]
            x0_old = x_arr[-2][0]
            x1_old = x_arr[-2][1]            
            ratio = 1 - f(x0_new,x1_new)/f(x0_old,x1_old)
        if n ==3:  
            x0_new = x_arr[-1][0]
            x1_new = x_arr[-1][1]
            x2_new = x_arr[-1][2]
            x0_old = x_arr[-2][0]
            x1_old = x_arr[-2][1]
            x2_old = x_arr[-2][2]            
            ratio = 1 - f(x0_new,x1_new,x2_new)/f(x0_old,x1_old, x2_old)
            
#        print(g_arr[-1])
#        print(x_arr[-1])
        if 0<ratio<1e-8  :
            finished = True
    return  x_arr[1], [x_arr]
    
def Grad2(f,x): # Finite difference 2d 
    eps = 1e-5
    grad = np.zeros(2)
    grad[0] = 0.5 *( f(x[0] + eps, x[1]) - f(x[0] - eps, x[1]) )/eps
    grad[1] = 0.5 *( f(x[0], x[1] + eps) - f(x[0], x[1] - eps) )/eps    
    return np.array(grad)

def Grad3(f,x): # Finite difference 3d 
    eps = 1e-5
    grad = np.zeros(3)
    grad[0] = 0.5 *( f(x[0] + eps, x[1], x[2]) - f(x[0] - eps, x[1], x[2]) )/eps
    grad[1] = 0.5 *( f(x[0], x[1] + eps, x[2]) - f(x[0], x[1] - eps, x[2]) )/eps
    grad[2] = 0.5 *( f(x[0], x[1], x[2] + eps) - f(x[0], x[1], x[2] - eps) )/eps  
    return np.array(grad)
    
def Gnext(gn, delta, gamma):
    a = (np.outer(delta,delta))/(np.dot(gamma,delta))
    b = (np.matmul( gn, np.matmul(np.outer(delta,delta),gn) ) )
    c =  np.matmul(np.matmul(gamma,gn),gamma)   
    return gn + a - b/c

#===================================

def UncDeriv3(t_min,m_min,cs_min,h=1e-7): # Assumes no mixing 
    deriv3_t = (NLL3(t_min + h, m_min, cs_min) - 2*NLL3(t_min, m_min, cs_min) + NLL3(t_min - h, m_min, cs_min) ) / h**2
    deriv3_m = (NLL3(t_min, m_min + h, cs_min) - 2*NLL3(t_min, m_min, cs_min) + NLL3(t_min, m_min - h, cs_min) ) / h**2
    deriv3_cs = (NLL3(t_min, m_min, cs_min + h) - 2*NLL3(t_min, m_min, cs_min) + NLL3(t_min, m_min, cs_min - h) ) / h**2
    return np.sqrt(1 / deriv3_t), np.sqrt(1 / deriv3_m), np.sqrt(1/deriv3_cs)


def test2(x,y):
    return x**2 + y**2 

def test3(x,y,z):
    return x**2 + y**2 + x**2 + x*y
