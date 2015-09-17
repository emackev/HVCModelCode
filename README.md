# HVCModelCode
Code for model in Okubo et al.

Emily Mackevicius 7/18/2015, based on Hannah Payne's code which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo Okubo. 

## To run 1 iteration of the model:
###HVCIter
See this file for step-by-step model dynamics and learning.

## To get paper figures:
### AlternatingDifferentiation
Code to generate figure 5 a-f, which shows alternating seed neuron differentiation, from subsong through protosyllable stage through splitting.

Relies on: 
####HVCIter
Runs 1 iteration of the model.  See this file for step-by-step model dynamics and learning.

####plotSubsong
Plotting function for subsong network diagram and raster.

####plotHVCnet
Plots network diagram.

####plotAlternating
Plots network activity

####findLatency
Called by plotting functions, tests what neurons participate in each syllable, and at what latencies.

###BoutOnsetDifferentiation
Code to generate Extended Data Figure a-d, which shows bout onset differentiation.

Relies on: 
####HVCIter
Runs 1 iteration of the model.  See this file for step-by-step model dynamics and learning.

####plotHVCnet_boutOnset
Plotting function for bout onset network diagram and raster.
  
####findLatency
Called by plotting functions, tests what neurons participate in each syllable, and at what latencies.

###SigLatDistOverDev.m
Latency distribution over development

###RunHVC_boutOnsetElement.m
Bout onset element differentiation (EDF e-h)

###RunHVC_split_intoThree.m
Motif learning, (EDF i-k)

###runHVC_split_movie.m
Alternating differentiation movie for supp

###boutOnsetDifferentiation_movie.m
Bout onset movie for supp

###RunHVC_boutOnset_Movies
Bout onset movie to visualize (sorted) weight matrix and activity over development, and multi-dimensional scaling of weight matrix
