# Complexity Project 


Author:  
Suhail Shoaib Mall  
Physics Department, Imperial College London  
CID 01392009  


Required non-standard packages:
* tqdm
* pickle 

Python files to run (in code folder):
* oslo_model.py
* heights.py
* avalanches.py
* get_pickles.py

## Example Usage 

### __OsloSim object:__ ### 
* Simulations can be run manually through the oslo_range.py file:

> s = OsloSim()

> s.Trun(100000)

* Will initialise and run a simulation for 100000 timesteps


### __heights.py__ ###

* Functions 'FindSteady' and 'SmootH' were used to generate simulations and data, but these have been stored for ease of use.
* Use GetParams() and GetSims() functions to retrieve data and simulations 

* Typically run functions in order that they are written with the output of one function as the input of the next
* Example of run to generate all data given in file :  
> DistributionQuestions()


### __avalanches.py__ ###

* Similar to heights.py in that data has already been run 
* Run functions in order of writing with out put of one as input of next
* Can run example function at end: 
> AvalancheExample() 


### __get_pickles.py__ ### 

* Use to load pickled simulations that have been run until they have reached steady state. 
* Also have pickled height and avalanche data 
