# Computational Physics Project 


Author:  
Suhail Shoaib Mall  
Physics Department, Imperial College London  
CID 01392009  


### __BACKGROUND:__ ###

* This project is to find the best-fit parameters to model the survival probability of a muon neutrino undergoing oscillation between neutrino flavours. 
* The total observed muons is then the expected number of observations multiplied by the probability of muon to muon oscillation. 
* This probability is dependent on the neutrino energy, as well as several mixing parameters.
* Accounting for this probaility in the expected number of observations, the model can be fit to real binned energy data by maximising the likelihood function. 
* This function is the likelihood of the model correctly describing the real binned energy data. Equivalently, one can minimise the negative log-likelihood (NLL). 

#### __Function Minimisers__ ####
This project explored several methods of minimising functions, including:
* Univariate Minimisation 
* Gradient Descent 
* Quasi-Newton Minimisation 

#### __Model Parameters__ ####
The model was first fit using two parameters:
* \theta_{23} -- the mixing angle 
* \Delta m^2_{23} -- the mass-square difference,
then with a third:
* cs -- a constant of proportionality to linearly mulitply the oscillation probability with neutrino energy 

The 3-parameter model was tested using Gradient Descemt and Quasi-Newton minimsation, and the two results compared.


### PROJECT STRUCTURE: ###


#### __nll_calcs.py__ ####

* First cell reads in the data hosted by Imperial College and uses the author's credentials to access it. 

* Data was read in and saved as .txt files , so use second cell to read data in using np.loadtxt(). Data is given as: 

	* lambda_arr : the expected number of observations without neutrino mixing, binned according to neutrino energy. Multiply this by the nuetrino mixing probaility to get the expeccted number of observations with mixing. 
	* k_arr : the actual number of observations per energy bin. 

* Contains functions to:
	* calculte the oscillation probability, total expected observations (with mixing), and NLL for two and three parameters.
	* Hessian and covariance matrices, given a function and position. This is useful for finding the standard error in a variable assuming the NLL is close to Gaussian at the minimm point.
	* Various plotting functions


#### __parabolic_approximation.py__ ####

* Contains functions to minimise the function varying a single parameter by fitting a Lagrange parabola to 3 points and finding the minimum. This is repeated until the minimum along this axis is found. 
* The univariate minimiser iterates through and minimises along the different parameters, repeating until a local minimum is found. 
* There are several functions to find the uncertainty by modelling as a Gaussian, aprabola, and finding the curvature using finite difference methods. 


#### __gradient_descent.py__ ####

* Contains methods to minimise the function by making steps against the function's gradient at that point.
* This is compeltely general, for any function of any number of variables, not just this particular NLL.


#### __quasi_newton.py__ ####

* Contains methods to minimise the function by taking steps against the gradient, but this time also considering the curvature of the function 
* This is also compeltely general, for any function of any number of variables, not just this particular NLL.



### USAGE EXAMPLES 

####  __plots:__ ####
> PlotOscProbE()  -- no output  
> Hists()  -- no output  
> HeatMap()  -- outputs (fig,ax) so that user can draw minimiser steps on top of heatmap  

#### __Parabolic Approximation:__ ####
> [theta_min, m_min], [theta_tests, m_tests] = Univariate(guess_theta = [0.8,0.9,0.95],guess_m=[2e-3,2.5e-3,3e-3])  

#### __Gradient Descent:__ ####
> xmin_2d, xsteps_2d= GradDesc(f=NLL_2d, x=[0.8,3e-3], alpha=[1e-4,1e-9], eps=2*[1e-5])  
> xmin_3d, xsteps_3d= GradDesc(f=NLL_3d, x=[0.86,2.88e-3,1.], alpha=[1e-4,1e-9, 1e-3], eps=3*[1e-5])  

#### __Quasi_Newton:__ ####
> x_QN, G_QN, steps_QN = QN(f=NLL_3d, x=[0.86,2.88e-3,1.], alpha=[1e-4,1e-9, 1e-3], eps=3*[1e-5])  


### PROJECT QUESTIONS 

* Also have functions at the end of each file that that run the functions with the parameters used to generate the results in the project report. 
* quasi_newton.py also has a function to compare the results for the gradient descent and QN minimisers.

##### __parabolic_approximation.py:__ #####
* *univariate_2d()* -- runs the univariate minimiser 

##### __gradient_descent.py:__ #####
* *grad_min_2d()* -- runs gradient descent minimiser for the 2d NLL  
* *grad_min_3d()* -- runs gradient descent minimiser for the 3d NLL  

##### __quasi_newton.py:__ #####
* *qn_2d()* -- runs quasi-Newton minimiser for 2d NLL. Also returns the final G matrix and calculates the Hessian so that these two can be compared. 
* *qn_3d()* -- runs quasi-Newton minimiser for 3d NLL. Also returns the final G matrix and calculates the Hessian so that these two can be compared. 
* *min_3d()* -- runs gradient descent and QN minimisers and outputs the results to be compared.

