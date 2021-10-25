# Radiations-IntracladePredatorEvolution
This repository contains the MATLAB code associated to the paper entitled "The enrichment paradox in adaptive radiations: emergence of predators hinders diversification in resource rich environments".

To make the figures 1, 2A, and 3 please run the rutine MakeFigures123. This may take 5 minutes or more to run, depending on the machine.

This repository also includes the code to simulate the radiation using an individual-based model (see supporting information of the paper), and to read its outcome. The individual-based model implementation can be found in the folder IBM. In this folder, you can find 2 folders: runSimulation and VisualizeResults. The first contains the rutine to simulate the radiation. The simulation can be performed running the rutine IBMdiversification. To read the population statistics, you can run the rutine VisualizeResults in the folder with the same name. Currently, the rutine VisualizeResults reads the file 'Predator.txt' which is the outcome observed in figure S5B of the paper.
