networks
========

Generative models and ABC-SMC framework for biological networks

Files:

- orca.so/libcorca.so: Cython/C++ libraries to count orbits (method from Hočevar, Demšar. Bioinformatics. 2014;30(4):559-65.)
- modelling_functions.py: Graph generative models and other utilities
- loadFunctions.py: Functions necessary for performing ABC-SMC
- abcsmc.py: script to submit ABC-SMC jobs

Directories:

- orca: C++ and Cython code for graphlet counting. I would recommend doing as much as possible to avoid recompiling... It can be frustrating
- plotting: Utilities for plotting posterior densities as approximised by the algorithm

GUIDELINES:
===========

- If you want to perform ABC-SMC on networks:

You must have a generative model as a python function that returns a networkX graph object.
First edit the abcsmc_graphs function in loadFunctions.py to call your function.
Then you can simply run abcsmc.py with a the required and optional arguments (python abcsmc.py to see them).

Liepe et al., Bioinformatics. 2010 Jul 15;26(14):1797-9. For details on the meaning of different options.
Toni et al., J R Soc Interface. 2009 Feb 6;6(31):187-202. To see how the algorithm works.


- If you want to run the generative DD model or versions of it:

All required functions are in modelling_functions.py



- If you want to plot ABC-SMC results:

Run pickle2table.py on your pickle output from abcsmc.py

Edit the R script to include your file names, parameter names and ranges.