# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:27:26 2022

@author: suhai
"""
from nll_calcs import * 

#%%    1D MINIMISATION 
    
def MinPoint(x0,x1,x2, y0,y1,y2): 
    '''
    > Given three (x_i,y_i) data points, fit to a quadratic using Lagrange
    polynomials and return the x-value that gives the minimum
    -------
    PARAMETERS: 
    xi : three x values 
    yi : three y values s.t. the data point is given by y_i = f(x_i) 
    -------
    RETURNS 
    x_min : fitting a quadrativ to the three data points, this is the x-value 
    for which the minimum occurs 
    '''        
    return 0.5 * ( (x2**2 - x1**2)*y0 + (x0**2 - x2**2)*y1 + (x1**2 - x0**2)*y2) / ((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2)


def ParaMinTheta(guess_points = [0.8,0.9,0.95], m=2.4e-3):
    '''
    > Parabolic minimisation of NLL varying just the theta parameter 
    > Approximate as parabola and find minimum point, use lowest three points 
    to fit the next parabola. Repeat until finishing condition met 
    -------
    PARAMS
    guess_points : list of initial values to start the search for the minimum
    theta value that minimises NLL
    m : m_squared value to be held constant 
    -------
    RETURNS 
    xmin : the value of theta that minimises nll 
    [xi] : list of theta values that fit the final parabola -- used in 
    univariate minimiser as the first guess for the next minimisation step
    '''
    change = 1. # initialised so that it doesn't terminate immediately 
    x0,x1,x2 = np.sort(guess_points)
    xmin_old = x0
    
    while change > 1e-3:           
        xmin =  MinPoint(x0,x1,x2,NLL_2d(x0,m),NLL_2d(x1,m),NLL_2d(x2,m))

        points = [x0,x1,x2,xmin]
        nll_points = [NLL_2d(p,m) for p in points]
        
        for point in points:
            if np.max(nll_points) == NLL_2d(point,m):
                points.remove(point)
                break 
        x0,x1,x2 = np.sort(points) 
        change = np.abs(xmin_old - xmin)/xmin_old
        xmin_old = xmin
    return xmin,[x0,x1,x2]  # x_min doesn't necessarily belong to [x0,x1,x2]

def ParaMinM(guess_points = [2e-2,2.5e-3,3e-3], theta = 0.95):
    '''
    > Parabolic minimisation of NLL varying just the m parameter 
    > Approximate as parabola and find minimum point, use lowest three points 
    to fit the next parabola. Repeat until finishing condition met 
    -------
    PARAMS
    guess_points : list of initial values to start the search for the minimum
    m value that minimises NLL
    theta : theta value to be held constant 
    -------
    RETURNS 
    xmin : the value of m that minimises nll 
    [xi] : list of m values that fit the final parabola -- used in 
    univariate minimiser as the first guess for the next minimisation step
    '''
    ratio = 1.
    x0,x1,x2 = np.sort(guess_points) # make sure in order 
    xmin_old = x0
    
    while ratio > 1e-3: # Arbitrary hyper-parameter            
        xmin =  MinPoint(x0,x1,x2,NLL_2d(theta,x0),NLL_2d(theta,x1),NLL_2d(theta,x2))

        points = [x0,x1,x2,xmin]
        nll_points = [NLL_2d(theta,p) for p in points]
        
        for point in points:
            if np.max(nll_points) == NLL_2d(theta,point):
                points.remove(point)
                break    
        x0,x1,x2 = np.sort(points) 
        ratio = np.abs(xmin_old - xmin)/xmin_old
        xmin_old = xmin
    return xmin,[x0,x1,x2]  

#%%        Uncertainty in theta_{23}

def UncGauss(theta_min, m=2.4e-3, sig_guess=0.05):
    '''
    > Treat likelihood as Gaussian, then \sigma is approximately found as the 
    value where the log-likelihood decreases by an absolute value of 0.5
    > Seek and return the values either side of the minimum point that satisfy
    this using scipy.optimize.fsolve
    -------
    PARAMS 
    theta_min : theta value that minimises nll 
    sig_guess : First guess of sigma to start search using op.solve 
    -------
    RETURNS 
    sig_minus : std for theta < theta_min
    sig_plus : std for theta > theta_min
    '''
    nll_min = NLL_2d(theta_min, m)
    theta_minus = op.fsolve(lambda theta : NLL_2d(theta,m)-nll_min-0.5 , 
              theta_min-sig_guess, full_output=True)
    theta_plus  = op.fsolve(lambda theta : NLL_2d(theta,m)-nll_min-0.5 , 
              theta_min+sig_guess, full_output=True)
    
    sig_minus = theta_min - theta_minus
    sig_plus  = theta_plus - theta_min
    return sig_minus, sig_plus 


def UncDeriv(theta_min, m = 2.4e-3, h=1e-7):
    '''
    > Approximate the likelihood as a Gaussian and \sigma can be found through 
    the curvature, i.e. the second derviatve. 
    -------
    PARAMS
    theta_min : theta value that minimises nll
    h : used to calculate the derivative using the central finite 
    difference scheme 
    -------
    RETURNS 
    sigma : std assuming nll approximates a Gaussian 
    '''
    
    deriv2_t = (NLL_2d(theta_min + h, m) - 2*NLL_2d(theta_min, m) + NLL_2d(theta_min - h, m) ) / h**2
    return np.sqrt(1 / deriv2_t)

def UncPara(t_min,m,t0,t1,t2):
    '''
    > Assuming the minimum is parabolic, the final three points used to find 
    the minimum can be used to calculate the curvature of the function at the 
    minimum. This is then related to the uncertainty in the minimum point 
    -------
    PARAMS 
    t_min : theta value that minimises nll 
    m : m_squared value 
    t_i : the three values that approximate the nll as a parabola at its 
    minimum point 
    '''
    y0 = NLL_2d(t0,m)
    y1 = NLL_2d(t1,m)
    y2 = NLL_2d(t2,m)    
    return (2*(y0/((t0-t1)*(t0-t2)) + y1/((t1-t0)*(t1-t2)) + y2/((t2-t0)*(t2-t1))))**-0.5

#%%    UNIVARIATE MINIMISER 

def Univariate(guess_theta = [0.8,0.9,0.95],guess_m=[2e-3,2.5e-3,3e-3]):
    '''
    > Iterating through parameteres, hold the other parameters constant 
    and vary this parameter to find the value that minimises the NLL 
    > Iterate through parameters and repeat until finishing condition met
    > Plot NLL to inform starting guess values 
    -------
    PARAMS
    > guess_theta : list of values of theta to start the search from 
    > guess_m     : list of values of m to start the search from 
    RETURNS 
    [t_min, m_min] : theta and m values that minimise the NLL, pressented in a 
    list  
    [t_tests, m_tests] : list of theta and m values that minimise NLL at each 
    step, such that plotting these would show the steps that the minimiser took 
    '''
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
        
        nll = NLL_2d(t_min,m_min)
        ratio = 1 - nll/nll_old
        nll_old = nll      
        if  np.abs(ratio) < 5e-8:
            finished = True 
    return [t_min, m_min], [t_tests, m_tests]


#%%        PROJECT QUESTION  

def univariate_2d():
    x_min, x_arr = Univariate()
    print('NLL minimised in 2d for values theta=%.4e and m=%.4e'%(xmin[0],xmin[1]))


