# General
This repository contains a separate data folder with the data needed to run the scripts for figure 4, figure 5 and figure 6. In addition, there is a folder that contains the figures that the script should give as output. Each script contains the code for an individual figure.

# Fig4.m
The first section of the script loads all the required data from the data folder. Spike_mat contains the spike times of the sorted units. Burst_edges contains the start and end frames of the population bursts. Pop_act contains the population activity data and pop_act_times contains the corresponding frames. See materials and methods section for more information. Only data for L1 is included so the other entries in the ORGANOIDS variable remain empty. 

The second section of the script can be used to control the function that creates figure 4. The variables in capital are adjustable. The settings in the script are the settings used for the manuscript.

# Fig5.m
This script requires the following toolbox:
http://strategic.mit.edu/downloads.php?page=matlab_networks 

The first section of the script loads all the required data and allows the user to set the settings for the STTC analysis. Spike_mat contains the spike times of the sorted units. Coordinates contains the spatial coordinates of the sorted units on the array. spike_times contains a cell array with the spike times of the sorted units per cell. See materials and methods section for more information on the analysis parameters. 

The second section of the script processes the data and computes the STTC between all pairs of units and filters the pairs based on their spike time latency distributions. 
The third section of the script saves the computed values so that they can be used for figure 6 as well.
The fourth section of the script can be used to control the function that creates figure 5 panels A-C. The variables in capital are adjustable. The settings in the script are the settings used for the manuscript.
The fifth section of the script can be used to control the function that creates figure 5 panels D-F. The variables in capital are adjustable. The settings in the script are the settings used for the manuscript.

# Fig6.m
The first section of the script loads all the required data and allows the user to set the settings for the STTC analysis. sttc_whole contains the connectivity matrix computed in Fig5.m. Coordinates contains the spatial coordinates of the sorted units on the array. Only data for L1 is included so the other entries in the ORGANOIDS variable remain empty. 
The second section of the script processes the data. The largest component of the connectivity network is selected and all the network edges are extracted. See materials and methods section for more information.

The third section of the script can be used to control the function that creates figure 6. The variables in capital are adjustable. The settings in the script are the settings used for the manuscript.
