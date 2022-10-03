# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:57:53 2019

@author: suhailMall
"""
import numpy as np

#%%
def epsilon(prec = '64'): # choose precision
    '''
    > Function to find the machine accuracy (epsilon) for numbers stored using
    the IEEE-754 standard for floating point representation
    > Machine accuracy is the smallest numebr than can be added to 1 without 
    loss of accuracy
    > It can be calculated as 2^{-n}, where n is the number of bits in the 
    mantissa   
    
    PARAMETERS
        > prec : precision --  the numpy datatype to be tested, either 32-bit, 
        64-bit, or extended precision. 
        Possible values: '32', '64', 'ext' respectively 
        
    RETURNS 
        > Denary value equal to the machine accuracy 
    '''
    if prec == '32': # make sure all numbers are of same precision 
        a0 = np.float32(1)
        mult = np.float32(0.5)
    elif prec == '64':
        a0 = np.float64(1)
        mult = np.float64(0.5)
    elif prec =='ext':
        a0 = np.longdouble(1)
        mult = np.longdouble(0.5)
    
    eps = a0*mult # initialise eps
    
    while a0+eps != a0:
        eps*=mult
    return eps/mult # gone one step too far --> divide by mult 


def small_rep(prec='64'): # smallest representable number 
    '''
    > Function to find the smallest represntable number for an un-normalised 
    IEEE-754 float.
    > The smallest representable number is 2^{-(n+b)} where n is the number 
    of bits in the mantissa and b is the bias.
    
    PARAMETERS
        > prec : precision --  the numpy datatype to be tested, either 32-bit, 
        64-bit, or extended precision. 
        Possible values: '32', '64', 'ext' respectively 
        
    RETURNS 
        > Denary value equal to the smallest representable number 
    '''
    if prec == '32':
        a0 = np.float32(1)
        mult = np.float32(0.5)
    elif prec == '64':
        a0 = np.float64(1)
        mult = np.float64(0.5)
    elif prec =='ext':
        a0 = np.longdouble(1)
        mult = np.longdouble(0.5)
        
    while a0  !=0.:
        a1=a0
        a0*=mult
    return a1

#%%

def q1(): # Returns results for all precisions 
    '''
    > No input parameters 
    > Returns the machine accuracy for single, double, and numpy extended 
    precision floats using the IEEE-754 standard. 
    '''
    single = epsilon('32')
    double = epsilon('64')
    ext = epsilon('ext')
    print('Single Precision: %.4e'%single)
    print('Double Precision: %.4e'%double)
    print('Extended Precision: %.4e'%ext)