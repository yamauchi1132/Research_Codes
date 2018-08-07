#############################
# check mechanical  energy  #
#############################

you can run this program as follow.

########################################################
python energy_check.py data1 data2 time.log (timestep) #
########################################################

Data1 is initial data such as sph_t0000.dat. Data2 is final data. 
In Data 2, the distance between two stars is the same as the distance in Data 1
Timestep, 3th argument, is the time data2 is made.
Data1, data2, and time.log is necessary arguments but timestep is option
If no timestep argument, timestep will be final time on the last row of the time.log automatically.
