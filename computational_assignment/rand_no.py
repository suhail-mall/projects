# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 17:56:38 2019

@author: suhailMall
"""
#import random as rand
import matplotlib.pyplot as plt
import numpy as np
from tqdm import trange
import tqdm
import random
import timeit

'''
This assignment was done in a notebook format as the question pertained to 
specific functions  
'''
#%%        Uniform deviate 
'''
> Validate the np.random.uniform distribution by considering its histogram 
'''
r = np.random.uniform(size=int(1e6)) # generate 10^6 uniform deviates
plt.figure()
plt.hist(r,bins=100, density = True, rwidth = 0.9, label = 'Histogram')
plt.hlines(1,0,1,label='PDF', color = 'red')
plt.xlabel('y')
plt.legend(loc='lower right')
#%%        TRANSFORMATION METHOD
# PDF: P(x) = 0.5cos(0.5x) , 0<x<pi
# Transform using y = F^{-1}(x), where F(x) = \int_{-\infty}^{x} P(z)dz
# then y is distributed according to P(x) 
# leads to transformation function y = 2*arcsin(x) 
# where x is given by uniform deviate 

# Use transformation method to convert uniform deviate to ~cos() dist. 
r = np.random.uniform(size=int(1e6))
y = [2 * np.arcsin(x) for x in tqdm.tqdm(r)] 
plt.figure()
plt.hist(y,bins = 100, density = True, rwidth = 0.9, label = 'Histogram')

# Test transformation method by comparing to ~cos() dist. 
a = np.arange(0.,np.pi,0.01) # x axis
b = [0.5*np.cos(0.5*x) for x in a] # PDF
plt.plot(a,b,'-', label = 'PDF')
plt.xlabel('y')
plt.legend()
#%%        REJECTION METHOD 
def PDF(x): 
    return (2/np.pi) * (np.cos(x/2))**2

def Comp(x): 
    return (2/np.pi) * np.cos(x/2)
#    return np.cos(x/2) 

def Comp_inv(x): # use to distribute according to comp funcn -- This is CDF
    return 2 * np.arcsin(x)

counter = 0
y=[] # returned array 
test_counter = 0 # used to count trials 

while counter <1e5:
    test_counter+=1 # increase each test 
    r = np.random.uniform() # uniform deviate from 0 to 1
    y_i = Comp_inv(r) # distributed according to comp funcn from 0-->pi
    p_i = np.random.uniform(0.,Comp(y_i))
    if PDF(y_i)>p_i:  
        y.append(y_i)
        counter+=1 
print(f'No. tests: {test_counter} ')

#plotting 
plt.figure()
a = np.arange(0.,np.pi,0.01)
b = [PDF(d) for d in a]
c = [Comp(d) for d in a]
plt.plot(a,b,'-', color = 'r', label = 'PDF')
plt.plot(a,c,'-', color = 'g', label = 'Comparison Function')
plt.hist(y,bins = 100, rwidth=0.9, density = True, label = 'Historgam')
plt.xlabel('y')
plt.legend()

#%%        Compare timing 

def time_rej(): # Just wrap code in a function to pass to timeit
    counter = 0
    y=[] # returned array 
    while counter <1e5:
        u = np.random.uniform() # uniform deviate from 0 to 1
        y_i = Comp_inv(u) # distributed according to comp funcn from 0-->pi
        p_i = np.random.uniform(0.,Comp(y_i))
        if PDF(y_i)>p_i:
            y.append(y_i)
            counter+=1
        return y

def time_trans():
    r = np.random.uniform(size=int(1e5)) 
    y = [2 * np.arcsin(x) for x in r] 
    return y 

print('Transformation method : %.4e s'%timeit.timeit(time_trans, number = 10))
print('Rejection method : %.4e s'%timeit.timeit(time_rej, number = 10))        
