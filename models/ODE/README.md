ODETestFull.py defines the ODE Model with noisy observations and adds an idle time to each model evaluation for a given runtime variance.
ODEWLogfiles.py does basically the same thing and saves all outputlogs. This is the newest version.
Parameter Inference is then performed repeatedly for different population sizes and for each of the different type of samplers.

ODEBoxPlotCreation.py evaluates the quality of the results generated in the Test
The .ipynb file can compute and visualize parallel efficiencies after running the Test several times on servers with different amounts of active workers.
