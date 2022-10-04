# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 16:57:42 2019

@author: suhai
"""
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from tqdm import tqdm


#%%        Check Matrix Dimension
def check(a):
    '''
    > Check that the input matrix is 2d and square (NxN)
    > If square, return N 
    > If not, raise Exception 
    
    INPUT 
        a : matrix A to be checked 
    RETURNS
        N : length of matrix of shape NxN 
        raise Exception : if A not square 
    ''' 
    a = np.array(a)
    shape = a.shape 
    if len(shape) !=2 or shape[0]!=shape[1]:
        #return None
        raise Exception ('Matrix is not square')
    return shape[0]

#%%        LU Decomp

def decomp(a):
    '''
    > Decompose matrix A into matrices L and U such that L.U = A
    > L is the lower matrix such that L_{ij}=0 for i<=j
    > U is the upper matrix such that U_{ij}=0 for i>j
    PARAMETERS 
        a : matrix to decompose 
    RETURNS 
        lu : decomposed matrix 
    '''
    n = check(a) 
    lu = np.zeros_like(a) #lu matrix is same dimension as a
    for i in range(n):
        for j in range(n): #count through elements
            if i >= j: #lower matrix 
                lu[i,j] = a[i,j] - np.sum([lu[i,k]*lu[k,j] for k in range(j)])
            elif i<j: #upper matrix 
                lu[i,j] = (a[i,j] - np.sum([lu[i,k]*lu[k,j] for k in range(i)])) / lu[i,i]
    return lu  
   

def l_extract(lu):
    '''
    > Given a LU matrix, return the lower part of the matrix 
    (including the diagonal) 
    PARAMETERS
        lu : LU matrix 
    RETURNS
        L : lower matrix 
    '''
    l = np.zeros_like(lu)
    n=len(lu)
    for i in range(n):
        for j in range(i+1):
            l[i,j]=lu[i,j]
    return l

def u_extract(lu):
    '''
    > Given a LU matrix, return the lower part of the matrix 
    > diagonal is set to 1. by convention 
    INPUT 
        lu : LU matrix 
    RETURNS
        u : upper matrix 
    '''
    u = np.zeros_like(lu)
    n = len(lu)
    for i in range(n):
        u[i,i]=1.
        for j in range(i+1,n):
            u[i,j] = lu[i,j]
    return u

#%%        Forward and Backward Substitutions 

def f_sub(l,b): 
    '''
    > Forward substitution to solve the matrix equation L.y = b 
    > L is a given lower matrix, b is a given vector 
    > y is the vector being solved for 
    > Can solve for y_i in series due to the triangular nature of L 
    PARAMETERS 
        l : lower matrix (shape NxN)
        b : vector (length N) 
    RETURNS
        y : vector that solves L.y = b 
    '''
    n=check(l) # Could be a repeat calculation here 
    if len(b) != n:
        raise Exception( 'Vector is Not Same Dimension as Matrix')
    y=np.zeros(n)
    for i in range(n):
        sum_term = 0.
        for j in range(i):
            sum_term += l[i,j]*y[j]
        y[i] = (b[i] - sum_term)/l[i,i]
    return y    

def b_sub(u,y): 
    '''
    > Backward substitution to solve the matrix equation U.x = y 
    > U is a given upper matrix, y is a given vector 
    > x is the vector being solved for 
    > Can solve for x_i in series due to the triangular nature of U 
    PARAMETERS 
        u : upper matrix (shape NxN)
        y : vector (length N) 
    RETURNS
        x : vector that solves U.x = y 
    '''
    
    n = check(u)
    if len(y) != n:
        raise Exception('Vector is Not Same Dimension as Matrix')
    x=np.zeros(n)
    for i in range(n-1,-1,-1): #count backwards
        sum_term = 0.
        for j in range(i+1,n):
            sum_term+=u[i,j]*x[j]
        x[i] = (y[i] - sum_term)
    return x

#%%        Solve Ax=b

def crout_axb(a,b):
    '''
    > A function using Crout's method to solve the matrix equation Ax=b, 
    where A is a given matrix, b is a given vector, and x is being sovled for 
    > Use LU decomposition on A so that L.(Ux)=b
    > Forward substition to solve L.y=b for y 
    > Backward substitution to U.x=y for x  
    PARAMETERS 
        a : matrix A (shape NxN)
        b : vector b (length N)
    
    RETURNS
        x : vector (length N)
    '''
    lu = decomp(a)
    l  = l_extract(lu)
    u  = u_extract(lu)
    return b_sub(u, f_sub(l, b))

#%%        Det and Inv

def det(a):
    '''
    > A function to calculate the determinant of a matrix using LU Decomp
    > Can show that det(A) is the product of diagonal elements of LU 
    PARAMETERS 
        a : matrix A to find det(A)
    OUTPUT 
        det(A)
    '''
    lu=decomp(a)
    return np.prod([lu[i,i] for i in range(len(lu))])

def inverse(a):
    '''
    > A function to find the inverse of matrix A using Crout's method 
    
    PARAMETERS 
        a : matrix A 
    RETURNS
        b : matrix A^{-1} -- the inverse of A 
    '''
    ide = np.zeros_like(a) #identity
    inv = np.zeros_like(a)
    for i in range(len(ide)): #index (i,j)
        ide[i,i]=1.        
        rhs = np.transpose(ide[:,i][None,:]) # turn the row vector into column
        inv[:,i] = np.transpose(crout_axb(a,rhs)) # picks column i 
    return inv 

#%%        Project Question 

a=np.array([[  3.,   1.,   0.,   0.,   0.],
       [  3.,   9.,   4.,   0.,   0.],
       [  0.,   9.,  20.,  10.,   0.],
       [  0.,   0., -22.,  31., -25.],
       [  0.,   0.,   0., -55.,  61.]])
    
b = np.array([2.,5.,-4.,8.,9.])

def q2():
    x = crout_axb(a,b)
    det_a = det(a)
    inv_a = inverse(a)
    print('x=', x)
    print('--------------------------------------------------------')
    print('Det(a) = ', det_a)
    print('--------------------------------------------------------')
    print('Inverse(a) = ',inv_a)