README for SSM2617, 01392009 Complexity Project

Required non-standard packages:
tqdm
pickle 

Python files to run (in code folder):
oslo_range.py
heights.py
avalanches.py
get_pickles.py
--------------------------------------------------------------------
Simulations:
Simulations can be run manually through the oslo_range.py file
e.g. 
s = OsloSim()
s.Trun(100000)

will initialise and run a simulation for 100000 timesteps

-------------------------------------------------------------------
heights.py

Functions 'FindSteady' and 'SmootH' were used to generate simulations and data, but htese have been stored for ease of use.
Use GetParams() and GetSims() functions to retrieve data and simulations 

Typically run functions in order that htey are writtenm with the output of one function as the input of the next
Example of run nto generate all data given in file = 'DistributionQuestions()'

-------------------------------------------------------------------
avalanches.py

Similar to heights.py in that data has already been run 
Run functions in order of writing with out put of one as input of next
Can run example function at end = 'AvalancheExample()'

-------------------------------------------------------------------
get_pickles.py 

Use to load pickled simulations that have been run until they have reached steady state. Also have pickled height and avalanche data 
