# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:23:00 2019

@author: suhailMall
"""

import scipy.fftpack as fft
import numpy as np 
import matplotlib.pyplot as plt 

#%%        Supporting Funcns

def gauss(t): # generate g(t) with specific params for question 
    return 1/np.sqrt(2*np.pi) * np.exp(-0.5*(t)**2)
def heavi(t):
    if 3.<t<5.:
        return 4. 
    return 0. 
    
def gauss_heavi(lower,upper,N):   
    # g = Gaussian funcn
    # h = Heaviside funcn
    # padded with N/2 0's after funcn  
    g = [gauss(t) for t in np.linspace(lower,upper, int(N/2))] + [0. for i in range(int(N/2))]  # second term is for padding 
    h = [heavi(a) for a in np.linspace(lower,upper,int(N/2)) ]+ [0. for i in range(int(N/2))] 
    return g,h

def ft(func):# Fast Fourier Transform 
    return fft.fft(func, n = len(func))

def inv(f_g,f_h,phase): # Find the prduct of the transformed functions and phase, then inverse transform back to time domain 
    product = [phase*f_h*f_g] # transformed arrays are still length N. T & N are defined in parent ffuncn so multiply there (constant)
    return fft.ifft(product)

def Phase(w,t_max): # Calculate phase to centre the function at 0 instead of +tmax
    return np.exp(1j*2*np.pi*w*t_max) # Due to definition of zero position



#%%        Convolution 
def conv(t_max=10,m=11, funcns=gauss_heavi): # t_max defines symmetrical range about 0. Unpadded is -tmax-->tmax. Then pad 2*tmax 
    T = 4*t_max 
    N=int(2**m)
    g,h = funcns(-t_max,t_max,N) # Generate functions and pad 
        
    f_g = ft(g) # Fourier Transform 
    f_h = ft(h)
    
    w = np.fft.fftfreq(N, d = T/N) # find w for phase. Need to do here as T & N are defined here
    w_phase = Phase(w,t_max)
    gh = inv(f_g,f_h,w_phase) * (T/N)  # delta t factor here as constant  
    return gh[0] # returns (N,1) dimension array 

def q4():
    m = 11 # defined here for plotting purposes 
    N = int(2**m)
    t_max = 10. 
    
    #For plotting 
    x = np.linspace(-t_max,3*t_max,N) # x-axis
    g,h = gauss_heavi(-t_max,t_max,N)  
    gh = conv(t_max,m,gauss_heavi)
    plt.plot(x,g,'--',label = 'g', color = 'green')
    plt.plot(x,h,'--',label = 'h', color = 'red')
    plt.plot(x,gh,'-',label = 'g*h', color = 'blue')
    plt.legend()
    plt.xlabel('t')
    plt.grid()
    
    
    
    

