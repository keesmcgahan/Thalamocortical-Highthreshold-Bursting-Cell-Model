# Thalamocortical-Highthreshold-Bursting-Cell-Model

Here we have the variety of different codes required to run the simulations which generate the different figures throughout the main text and supplement of the paper "Specific intrinsic currents of an LGN thalamocortical cell model allow for awake firing dynamics necessary to explain retinal transfer properties".

The first, most important of these files is landt.ode. THis file is what is meant to be loaded into xppauto in order to produce the bifurcation code. It also is the cleanest representation of the parameters, model equations and ODES without any additional code for other simulations. 

tcGAfit.m is Matlab code to run a genetic algorithm which fits a thalamocortical cell model against the features which describe low threshold and high threshold bursting. In particular the model is fit to reporudece the different frequencies, and intraburst intervals, and minimum voltage requirements. These can be adjusted to also include changes in variance of ISI spikte times. Additionally weights and spific desired values can be adjusted. Future simulations are based on one possible parameter set.

lmhtINPUTS.m is the primary matlab code for this study. In here it enables the used to simulate the neuron to steady state for any giving underlying dynamic regime, and then input different excitatory retinal signals. The strength and timing of these incoming spiking signals can be adjusted to be poisson or gamma distributed, and have any desired type of mean spikes per second desired. This code simulateously will output the membrane potential plots as well as the accompnaying histogram which showcase the delay times between ret/lgn spike pairs, as well as showcasing the histogram which shows the relationship between retinal ISI and output of an LGN spike.

THe final two codes are for the last section of the main document and its accompnaying supplmenetal figures. THese are ExcitatoryEntrainment.m and InhibitionsEntrainment.m. These two files enable the user to simulate the nueron in a desired dynamic system regime in response to a desired type of rhythmic stimulation. 
