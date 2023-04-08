# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:02:12 2022

@author: suhai
"""

from gradient_descent import * # To get the derivative function 

def Gnext(G_old, delta, gamma):
    '''
    > Uses the DPL algorithm to approximate the inverse Hessian matrix 
    > Found iteratively, updating G at position x_{n+1} based on G at x_n
    and the function and its derivatives at x_n and x_{n+1}
    ----------
    PARAMETERS 
    G_old : the matrix G at position x_n in parameter space 
    delta : f(x_{n+1}) - f(x_n) 
    gamma : \grad f(x_{n+1}) - \grad f(x_n)
    -------
    RETURNS
    G : approximate inverse Hessian at position x_{n+1}
    '''
    a = (np.outer(delta,delta))/(np.dot(gamma,delta))
    b = (np.matmul( G_old, np.matmul(np.outer(delta,delta),G_old) ) )
    c =  np.matmul(np.matmul(gamma,G_old),gamma)   
    return G_old + a - b/c

def QN(f, x, alpha, eps):
    '''
    > Finds the minimum of any function of n params by taking steps against the 
    gradient and considering the curvature of the function at that point 
    > Make steps according to: x_{n+1} = x_n - \alpha G \nabla f(x_n), where 
    G is approximately the inverse Hessian matrix 
    > G is found iteratively using the DFP algorithm 
    > See project report for more information 
    -------
    PARAMS 
    f : the function to be minimised 
    x : an n-dim vector of the first guess parameters to start the search 
    from 
    alpha : an n-dim vector of the step-size for each parameter 
    eps : an n-dim vector for finite differences along each parameter
    -------
    RETURNS  
    x_min : an n-dim vector the values of the parameters that minimise f 
    G : the (n,n) G matrix evaluated at the minimising position. This can be 
    used to calculate the curvature and standard error at this position 
    x_arr : an (N,n) array of the N steps taken by the minimiser. Each element 
    is the n-dim vector of the position of the minimiser at that step 
    '''
    finished = False # Initialise as False to start steps
    n=len(x)
    x = np.array(x) # Parameters 
    alpha = np.array(alpha) # Step distance for each param
    x_arr = np.array([x]) # List of steps 
    
    G = np.identity(n) # Initialise G as identity matrix 
    der =  der_vector(f, x, eps) # Find derivative at start point 
    
    while not finished:
        # Calculate step based on current position 
        # step = np.diagonal(np.outer(alpha, np.matmul(G,der))) 
        Gder = np.matmul(G,der)
        step = np.array([a*Gd for (a,Gd) in zip(alpha, Gder)])
        
        x_new = x - step # Move to next point
        der_new = der_vector(f, x_new, eps) # Calculate der at new point 
        delta = x_new - x # Use these to calculate G at new point 
        gamma = der_new - der 
        G_new = Gnext(G, delta, gamma)
        
        # Update G and der as these are now the "old" values after the step 
        G = G_new 
        der = der_new 
        
        ratio = 1 - f(*x_new)/f(*x) # Check finishing condition 
        if 0 < ratio < 1e-8:
            finished = True 
        else:
            x = x_new # If not finished, new position becomes old 
            x_arr = np.append(x_arr,np.array([x]),axis=0) # To track steps 
    return x_new, G, x_arr

#%%        PROJECT QUESTION 

def qn_2d():
    eps = 2*[1e-5]
    alpha = [1e-4,1e-9]
    xmin_2d, G, x_arr = QN(NLL_2d, [0.86,2.88e-3], alpha , eps)
    print (f'NLL minimised in 2d for values theta=%.4e,  and m=%.4e '%(xmin_3d[0],xmin_3d[1],xmin_3d[2]))
    print('Final G Matrix:')
    print(G)
    H = Hessian(NLL_2d, xmin_2d, eps)
    print('Hessian Matrix at final point: ')
    print(H)



def qn_3d():
    eps = 3*[1e-5]
    alpha = [1e-4,1e-9, 1e-3]
    xmin_3d, G, x_arr = QN(NLL_3d, [0.86,2.88e-3,1.], alpha , eps)
    print (f'NLL minimised in 3d for values theta=%.4e, m=%.4e, and cs=%.4e '%(xmin_3d[0],xmin_3d[1],xmin_3d[2]))
    print('Final G Matrix:')
    print(G)
    H = Hessian(NLL_3d, xmin_3d, eps)
    print('Hessian Matrix at final point: ')
    print(H)

#%%        COMPARE GRADIENT MINIMISER AND QN 

def min_3d():
    alpha = [1e-4,1e-9, 1e-3]
    eps = 3*[1e-5]
    xgrad, steps_grad= GradMin(NLL_3d,[0.86,2.88e-3,1.], alpha,eps)
    x_QN, G_QN, steps_QN = QN(NLL_3d, [0.86,2.88e-3,1.], alpha, eps)
    
    print('NLL minimised by parameters:')
    print('Gradient Minimisation: ', xgrad)
    print('Quasi-Newton: ', xQN)
    
    
    
    