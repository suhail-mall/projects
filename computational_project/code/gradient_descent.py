# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:17:45 2022

@author: suhai
"""
from nll_calcs import * 
#%%        GRADIENT MINIMISATION 

def der_vector(f,x,eps):
    '''
    > Find the first derivative of a function at a given position 
    -------
    PARAMS 
    f : a function of any n-dim variables 
    x : an n-dim vector describing a position to evaluate f at 
    eps : an n-dim vector for finite differences along each parameter 
    -------
    RETURNS 
    der : an n-dim derivative vector of f at x 
    '''
    eps_mat = np.diag(eps) # To add eps[i] to the i'th param 
    n = len(x)
    eps_mat = np.diag(eps) # To add eps[i] to the i'th param 
    der = np.zeros_like(x)
    
    for i in range(n): # iterate through params to calc der along each param using central diff 
        d_xi_plus  = x + eps_mat[i] # Add eps[i] to the ith param for finite diff
        d_xi_minus = x - eps_mat[i] 
        der[i] = (f(*d_xi_plus) - f(*d_xi_minus))/(2*eps[i]) # calculate der
    return der 

def GradDesc(f, x, alpha, eps):
    '''
    > Finds the minimum of any function of n params by taking steps against the 
    gradient 
    > Make steps according to: x_{n+1} = x_n - \alpha \nabla f(x_n)
    > See project report for more information 
    -------
    PARAMS 
    f : the function to be minimised 
    x : an n-dim vector of the first guess variables to start the search 
    from 
    alpha : an n-dim vector of the step-size for each parameter 
    eps : an n-dim vector for finite differences along each parameter
    -------
    RETURNS  
    x_min : an n-dim vector the values of the variables that minimise f 
    x_arr : an (N,n) array of the N steps taken by the minimiser. Each element 
    is the n-dim vector of the position of the minimiser at that step 
    '''
    finished = False # Initialise as False to start steps 
    x = np.array(x) # variables 
    alpha = np.array(alpha) # Step distance for each param
    x_arr = np.array([x]) # List of steps 
    
    while not finished: 
        der = der_vector(f, x, eps)# Calculate der vector 
        step = np.diagonal(np.outer(alpha,der)) # step[i] = alpha[i]*der[i]
        x_new = x - step # make step 
        ratio = 1 - f(*x_new)/f(*x) # check finishing condn 
        if 0 < ratio < 1e-8:
            finished = True 
        else:
            x = x_new
            x_arr = np.append(x_arr,np.array([x]),axis=0)
    return x_new, x_arr

#%%        PROJECT QUESTION 

def grad_min_2d():
    xmin_2d, xsteps_2d= GradDesc(NLL_2d,[0.8,3e-3],[1e-4,1e-9],2*[1e-5])
    print (f'NLL minimised in 2d for values theta=%.4e and m=%.4e'%(xmin_2d[0],xmin_2d[1]))
    c = cov(Hessian(NLL_2d, xmin_2d, [1e-4,1e-4]))
    print ('Covariance Matrix : ')
    print(c)
        
def grad_min_3d():
    xmin_3d, xsteps_3d= GradDesc(NLL_3d,[0.86,2.88e-3,1.], [1e-4,1e-9, 1e-3],3*[1e-5])
    print ('NLL minimised in 3d for values theta=%.4e, m=%.4e, and cs=%.4e '%(xmin_3d[0],xmin_3d[1],xmin_3d[2]))
    c = cov(Hessian(NLL_3d, xmin_3d, 3*[1e-4]))
    print ('Covariance Matrix : ')
    print(c)
    