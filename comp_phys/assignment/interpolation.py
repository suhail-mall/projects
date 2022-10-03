# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 23:13:37 2019

@author: suhai
"""
import numpy as np
import matplotlib.pyplot as plt
from crout import *

#%%        Gausian function (used for testing)
def gauss(x,x0,sig):
    '''
    > Gaussian function for a single variable 
    PARAMETERS
        x   : variable 
        x0  : mean 
        sig : std 
    Returns
        G(x;x0,sig)
    '''
    return 1/(np.sqrt(2*np.pi)*sig) * np.exp ( - ((x-x0)**2) /(2*(sig)**2) )

#%%

def xf_check(x,f):
    '''
    > Check that data points (x, f(x)) are same length 
    > x and f given parameterically, i.e. data point i given by (x[i],f[i])
    ----------
    PARAMETERS 
        x : x data 
        f : f(x) data 
    RETURNS 
        n : if x and f same length, return number of data points 
        else, raise Exception   
    '''
    n = len(x)
    if n!= len(f):
        raise Exception('X and F are of different dimension')
    return n

def lin_interp(x_data,f_data,x): 
    '''
    > function to linear interpolate data (x, f(x))
    > Given some {x,f}_i data and a point x, return an interpolation of f(x)
    ----------
    PARAMETERS 
        x_data : x data array 
        f_data : f data array of same length as x 
    RETURNS 
        f : f(x) given by interpolating between two given data points, at a 
        given x point 
    '''
    n = xf_check(x_data,f_data)
    if not x_data[0] <= x <= x_data[-1]:
        raise Exception('Cannot Interpolate for point outside data range: %s -> %s'%(x_data[0],x_data[-1]))
    for i in range(n):
        if x_data[i] >= x : # Find the two given data points that x is between 
            m0 = i-1
            m1 = i
            break # Don't want to keep counting
    return  ( (x_data[m1] - x)*(f_data[m0]) + (x-x_data[m0])*(f_data[m1]))/(x_data[m1] - x_data[m0])


def lin_plot(x,f, plot ):
    '''
    > Linear interpolate the given data at a resolution 100x that of the 
    given data
    ----------
    PARAMETERS 
        x : x data
        f : f data 
    RETURNS 
        x_interp : the points that were sampled 
        f_lin : the function at the x_interp points using linear interpolation
    '''
    x_min = np.min(x)
    x_max = np.max(x)
    x_interp = np.linspace(x_min,x_max,len(x)*100)
    
    f_lin = [lin_interp(x,f,a) for a in x_interp]
    
    if plot:
        plt.figure()
        plt.plot(x,f,'x', color = 'red', mew=1.5, label = 'data')
        plt.plot(x_interp,f_lin, '-', color = 'orange', label = 'linear interpolation')
        plt.grid()
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('f')
    
    return x_interp, f_lin 
#%%        Set up the matrix eqn
def mat_eqn(x,f): 
    '''
    > Set up the matrix equation to solve for the second derivative 
    at each data point.
    > See project report for details on this. 
    > Using a natual spline so the second der. of the first and final data 
    points are set to 0. They are not being solved for so remove these from 
    the matrix eqn. 
    ----------
    PARAMETERS 
        x : x data array 
        f : f data array of same length as x 
    RETURNS 
        mat : matrix that acts on the derivative vector, der.
        rhs : a vector that depends on the data around a point x_i.
        mat and rhs defined so that   mat.der=rhs , where der is the vector of 
        second derivatives at points x_i
    '''
    n = len(x)
    m=n-2 # cut top and bottom rows off
    mat=np.zeros([m,m])
    rhs = np.zeros(m)
   
    for i in range(1,m-1): # Not using f"_0 and f"_n so first and last rows only consider 2 derivatives --> solve serparately
        j = i+1 #shift indices up by 1 - cut out the first row
        mat[i,i-1] = (x[j] - x[j-1]) / 6
        mat[i,i] = (x[j+1] - x[j-1]) / 3
        mat[i,i+1] = (x[j+1] - x[j]) / 6
       
        rhs[i] = (f[j+1] - f[j])/(x[j+1]-x[j]) - (f[j]-f[j-1])/(x[j]-x[j-1])
    
    # Need to manually set points neighbouring the end points due to natural 
    # spline 
    # First and last rows depend on x_0 and x_n    
    rhs[0] = (f[2]-f[1])/(x[2]-x[1]) - (f[1]-f[0])/(x[1]-x[0]) #i-1 = 0
    rhs[-1] = (f[-1] - f[-2])/(x[-1]-x[-2]) - (f[-2]-f[-3])/(x[-2]-x[-3])
   
    mat[0,0] = (x[2]-x[0])/3   # f"_1 depends on x_0 and x_2
    mat[0,1] = (x[2]-x[1])/6
    mat[-1,-2] = (x[-2]-x[-3]) / 6 #bottom right corner -- i = m-1 (+1)
    mat[-1,-1] = (x[-1] - x[-3]) / 3 # one to the left of it -- i = m-1 (+1)
    #(i+1) = x[-1], (i) = x[-2], (i-1) = x[-3]
    return mat, rhs

def deriv(x,f):
    '''
    > Generate and solve matrix equation for the second der. vector 
    > Then insert 0. in first and last positions for the natual spline 
    ----------
    PARAMETERS 
        x : x data array 
        f : f data array of same length as x 
    RETURNS 
        der : array of second derivatives at the data points     
    '''
    mat,rhs = mat_eqn(x,f) # generate matrix eqn
    der = crout_axb(mat,rhs) # solve for f"_1 --> f"_{n-1}
    der=np.insert(der,0,0.) # insert 0 to top and bottom for f"_0 and f"_n
    der=np.insert(der,len(der),0.)  
    return der # list of derivatives 


def spline_point(x,f,x_user,der): 
    '''
    > Given the data points and their second derivatives, use a cubic spline to 
    calculate the function at a user-given point 
    ----------   
    PARAMETERS 
        x      : 1d x data array 
        f      : 1d f data array of same length as x 
        x_user : the x point that the user wants to calcualte f at 
        der    :  the array of second derivatives at each x point 
    RETURNS 
       f_spline : the value of f(x_user) found using the spline   
    '''
    n = xf_check(x,f)
    for i in range(n):
        if x[i] >= x_user :
            m0 =  i-1
            m1 = i
            break # Don't want to keep counting    
    a = (x[m1] - x_user) / (x[m1] - x[m0]) # See project report for details 
    b = 1 - a                              # on deriving these relations
    c = (1/6) * (a**3 - a ) * (x[m1] - x[m0])**2          
    d = (1/6) * (b**3 - b ) * (x[m1] - x[m0])**2
    return a*f[m0] + b*f[m1] + c*der[m0] + d*der[m1]


def spline_plot(x,f, plot=True): 
    ''''
    > Use a cubic spline to interpolate the given data at a resolution 100x 
    that of the given data
    > Find the second derivative for all data points using Crout's method 
    then use this without recalculating 
    ----------
    PARAMETERS
        x : x data array 
        f : f data array of same length as x 
        plot : Boolean value of if to plot data and interpolant 
    RETURNS 
        x_interp : the x data array at higher resolution 
        y_spline : f(x_interp) array using a cubic spline 
    '''
    x_min = np.min(x)
    x_max = np.max(x)
    x_interp = np.linspace(x_min,x_max,len(x)*100)
    
    der = deriv(x,f)
    y_spline = [spline_point(x,f,a,der) for a in x_interp]
    
    if plot:
        plt.figure()
        plt.plot(x,f,'x', color = 'red', mew=1.5, label = 'data')
        plt.plot(x_interp,y_spline,'-', color = 'blue', label = 'Cubic Spline')
        plt.grid()
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('f')
    return x_interp, y_spline
 
    
 #%%
    
def q3():
    '''
    > Interpolate some given data at 100x resolution using a linear interpolant 
    and a cubix spline. 
    '''
    x=np.array([-2.1,-1.45,-1.3,-0.2,0.1,0.15,0.9,1.1,1.5,2.8,3.8])
    f=np.array([0.012155,0.122151,0.184520,0.960789,0.990050,0.977751,
                0.422383,0.298197,0.105399,3.93669e-4,5.355348e-7])

    # Run each point through interpolator for points in data range then plot 
    der = deriv(x,f) # Run this once at start then pull required derivatives 
    x_min = np.min(x)
    x_max = np.max(x)
    x_interp = np.linspace(x_min,x_max,len(x)*100)
    y_spline = [spline_point(x,f,a,der) for a in x_interp]
    y_lin = [lin_interp(x,f,a) for a in x_interp]
    plt.figure()
    plt.plot(x,f,'x', color = 'red', mew=1.5, label = 'data')
    plt.plot(x_interp,y_spline,'-', color = 'blue', label = 'Cubic Spline')
    plt.plot(x_interp,y_lin, '-', color = 'orange', label = 'linear interpolation')   
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('f')
    





