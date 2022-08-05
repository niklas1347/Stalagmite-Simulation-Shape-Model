# niklas1347-Stalagmite-Simulation-Shape-Model
The Shape Model allows the user to simulate the formation of stalagmites with 4 input parameters: calcium concentration of the water drop, drip interval, cave temperature, and cave carbon-dioxide. The output Parameters are growth rate, equilibrium radius, excessive calcium concentration and the respective depth of the values. The Shape model is used at the example of the So-1 in the Sofular Cave, Turkey in the paper which soon will be puplished. More background information can be found there. If the original data for this paper is needed please contact the corresponding author. Please cite this paper if this code is used for further studies.

All codes are available in jupyter notebook and phyton.

auxiliary_funcionts: additional functions which are needed to run the Shape Model
plot_referenz_and_modify_stalagmite: algortihm to plot a stalagmite through the FLOW-Equations, additional function which starts growing in an equilibrium and some variations of these functions
plot_comparsion_models: algorithm to plot the stalagmite with equal axises, additional functions to plot just the upper line for the FLOW-, GAUSS- and EXP-Equations
plot_changes_in_cave: algorithm plot the stalagmite in the cave and to detect changes, addtional functions to get the values for a changing CO2 concentration in the cave
