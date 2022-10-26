# Computational Physics Assignment 

Author: Suhail Shoaib Mall  
	Physics Department, Imperial College London  
	CID 01392009  



### __Question 1 - Machine Accuracy (float.py):__ ###

* The Machine Accuracy code is run by simply calling the q1() function. 
* This runs the epsilon(prec) function for all available precisions.
* The returned values are the machine epsilon values for datatypes:
	* np.float32
	* np.float64
	* np.longdouble

### __Question 2 - Matrix Methods (crout.py):__ ###

* Calling the function q2() is used to retrieve answers to question.
* However decomposition of any matrix a is done by simply calling decomp(a).
* Solving equation Ax=b for x is done by calling b_sub(A,b).
* Calling these outside the function q2() is useful as matrix lu and vector x can be viewed using the Spyder variable explorer.

* Determinant and Inverse of matrix A are found by calling det(a) and inverse(a) respectively. 


### __Question 3 - Interpolation (interpolation.py)__ ###

* Calling the function q3() plots of given data (x,f), the linear interpolation, and the cubic spline.
* For array data for linear iterpolator or cubic spline, call lin_plot(x,f,False) or call spline_plot(x,f,False) respectively.
* This returns x_interp array and interpolation array of same length.
* Third argument is whether or not to plot.


### __Question 4 - Fourier Transforms (fourier.py)__ ###

* Call q4() to plot functions g and h and convolution (g*h) for m=11 and t_max = 10. (2^10 points in range -10 ---> 10)
* For array of (g*h) for arbitrary (t_max,m) use gh = conv(t_max,m)
* This can be plotted against axis in t-domain:	
> x = np.linspace(-t_max,3*t_max,N)

> plt.plot(x,gh)


### __Question 5 - Random Numbers (rand_no.py)__ ###

* Run individual cells for distributions (cells are titled)
* The final cell is run to measure the time taken to run each generator 