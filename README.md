# Epi_SSA
# Epi_SSA
This repository contains MATLAB scripts for Stochastic Simulation Algorithms for Epidemics. It was built for my Master's Project *Differential Equations vs Stochastics Simulations Algorithms for Epidemics*, with the supervisor of Dugald Duncan. 

##Introduction
We use [Stochastic](https://en.wikipedia.org/wiki/Stochastic_simulation) [Simulation](https://en.wikipedia.org/wiki/Gillespie_algorithm) Algorithms to simulate outbreaks of Measles and Squirrel Pox. WWe also solve the equivalent SIR ODEs to compare the two techniques. This repository contains scripts for basic SIR outbreaks, to simulation of more complex models with heterogenous mixing and small work connections. 


##Basic SSA Algorithms
`IR_MODEL.M` simulates a group of individuals recovery from a disease. It simulates this stochastically, as well as solving an ODE model and compares the two outputs. 

`SIR_GILLESPIE.M` simualtes an outbreak of a disease in a large group of individual tracking the number of sucsceptibles, infected and recoverds. `SIR_ODE.m` is the accompanying function file containing the equivalent ODE model. 

We extend this model further in `Measeles_SIR.m` introducing birth and death processes. This script simulates periodic outbreak of mealses in a large population. It compares the ODE, which tends to a fixed point, againts the SSA, where the stochastic process means the outbreak can oscilate without tending to a fixed point. 

##Squirrel Pox Simulations
Squirrel Pox is a disease which is deadly to red squirrels, but not to greys.  These scripts are concerned with the introduction of squirrel pox into an island of red squirrels. In `SQUIRREL_POX_MODEL.m` and the accomponying fucntion `SQUIRREL_POX_ODE.m`, we conisider an island split up into a chain of boxes, with each box having a population of squirrels and a carrying capacity. We allow squirrels in each box to mix with their neighbours, which acts as the main method for disease transmission. We also allow suqirrels to disperse, meaning they can permemantly move to another randomly chosen box. 

The idea of dispersing to another box can be extened by implementing a small world network, where certain boxes are more strongly connected that others, giving a higher rate of dispersion. This is implemented in `SQUIRREL_POX_SW.m`, which also contains a scenrio where there is a barrier on the island which few squirrels can cross.  
