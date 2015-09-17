# HVCModelCode
Code for model in Okubo et al.

Emily Mackevicius 7/18/2015, based on Hannah Payne's code which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo Okubo. 

## To run 1 iteration of the model:
###[HVCIter](https://github.com/emackev/HVCModelCode/blob/master/HVCIter.m)
See this file for step-by-step model dynamics and learning.

## To get paper figures:
###[AlternatingDifferentiation](https://github.com/emackev/HVCModelCode/blob/master/AlternatingDifferentiation.m)
Code to generate figure 5 a-f, which shows alternating seed neuron differentiation, from subsong through protosyllable stage through splitting.

Relies on: 

[HVCIter](https://github.com/emackev/HVCModelCode/blob/master/HVCIter.m)
Runs 1 iteration of the model.  See this file for step-by-step model dynamics and learning.

[plotSubsong](https://github.com/emackev/HVCModelCode/blob/master/plotSubsong.m)
Plotting function for subsong network diagram and raster.

[plotHVCnet](https://github.com/emackev/HVCModelCode/blob/master/plotHVCnet.m)
Plots network diagram.

[plotAlternating](https://github.com/emackev/HVCModelCode/blob/master/plotAlternating.m)
Plots network activity

[findLatency](https://github.com/emackev/HVCModelCode/blob/master/findLatency.m)
Called by plotting functions, tests what neurons participate in each syllable, and at what latencies.

###[BoutOnsetDifferentiation](https://github.com/emackev/HVCModelCode/blob/master/BoutOnsetDifferentiation.m)
Code to generate Extended Data Figure a-d, which shows bout onset differentiation.

Relies on: 

[HVCIter](https://github.com/emackev/HVCModelCode/blob/master/HVCIter.m)
Runs 1 iteration of the model.  See this file for step-by-step model dynamics and learning.

[plotHVCnet_boutOnset](https://github.com/emackev/HVCModelCode/blob/master/plotHVCnet_boutOnset.m)
Plotting function for bout onset network diagram and raster.

[findLatency](https://github.com/emackev/HVCModelCode/blob/master/findLatency.m)
Called by plotting functions, tests what neurons participate in each syllable, and at what latencies.

###[HVCModelCode](https://github.com/emackev/HVCModelCode/blob/master/HVCModelCode.m)
Generates Figure 5 a-d and Extended Data Figure a-d. All in one .m file, including functions it depends on, for posting with the paper. 

###[SigLatDistOverDev.m](https://github.com/emackev/HVCModelCode/blob/master/SigLatDistOverDev.m)
Latency distribution over development (Fig 5 e).

###[RunHVC_boutOnsetElement.m](https://github.com/emackev/HVCModelCode/blob/master/RunHVC_boutOnsetElement.m)
Bout onset element differentiation (EDF e-h)

###[RunHVC_split_intoThree.m](https://github.com/emackev/HVCModelCode/blob/master/RunHVC_split_intoThree.m)
Motif learning, (EDF i-k)

###[runHVC_split_movie.m](https://github.com/emackev/HVCModelCode/blob/master/runHVC_split_movie.m)
Alternating differentiation movie for supp

###[boutOnsetDifferentiation_movie.m](https://github.com/emackev/HVCModelCode/blob/master/boutOnsetDifferentiation_movie.m)
Bout onset movie for supp

###[RunHVC_boutOnset_Movies](https://github.com/emackev/HVCModelCode/blob/master/RunHVC_boutOnset_Movies.m)
Bout onset movie to visualize (sorted) weight matrix and activity over development, and multi-dimensional scaling of weight matrix
