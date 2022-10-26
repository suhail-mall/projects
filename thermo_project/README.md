# Thermodynamics Project 



Author:  
Suhail Shoaib Mall  
Physics Department, Imperial College London  
CID 01392009  

This project is designed to simulate a non-ideal gas using hard sphere balls.



NOTE: 
* Must have following packages installed
	* _tqdm_
	* _numpy_ 
	* _scipy_
	* _scipy.optimise_
	* _matplotlib.pyplot_

* All of these packages are installed on the computers in Blackett Laboratory computer suite but if not, they can be installed easily through Anaconda or using pip install. If don't want to use TQDM, comment out import line and change trange to range.
      

## Project Structure:

* To take data use modules within TESTING folder.
* To view simulation animation use nBall module.

* Core modules:
	* initBall
	* nBall

### __Testing:__ ###
* Probes used to measure thermodynamic properties of a simulation
* Runs experiments by initalising simulations with different initial parameters 
	* distributionTest
	* presN
	* presTempContinuous
	* presTempReset
	* presTime
	* presV
	* thermoCalc

* Results:
	* Contains .png figures of plots of data
	* Data also included for most experiments
	* Naming convention inconsistent but semi-coherent
	* Each file contains array of just one variable (not both variables in each experiment)


## Use and Examples:

### __Initialising a Ball:__ ###
* Initialise ball with parameters:
        
	* m= mass of ball, float
        
        * r= radius of ball, float
        
        * pos = position of ball; list of length 2, each element of type float
        
	* vel = velocity of ball; list of length 2, each element of type float

* Can use collide() function as self.collide(other)
* Takes other Ball as input, e.g:

> b1=Ball(1.,0.1,[2.,0.],[-1.,0.])
> b2=Ball(1.,0.1,[-2.,0.],[11.,0.])
> b1.collide(b2)          


### __Initalising a Simulation:__ ###
* Initialise a Sim class using with a float = radius of Container

* Typically use Sim.generateBalls() to populate Container.
* Takes parameters:
	* n = number of Balls, int
	* m = mass of each Ball, float
	* r = radius of each Ball, float
	* v_bottom = lower bound of velocity-component
	* v_top = upperbound of velocity-component


* e.g. of initialisation and use:

> s=Sim(10.)
> s.generateBalls(100,1.,0.1,-10.,10.)
> s.Run(2000,True)


### __thermoCalc and distributionTest:__ ###

* Functions take Sim class and return thermodynamic quantity from system.

* Use example of presN module:

pN function takes parameters (cont,n,m,r,v_bottom,v_top), all defined above, e.g: 

> a=pN(10.,100,1.,0.1,-10.,10.)
> plt.plot(a[0],a[1],'x')






							
