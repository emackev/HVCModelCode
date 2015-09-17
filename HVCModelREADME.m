%% Code for HVC model in Okubo et al

% Emily Mackevicius 7/18/2015, based on Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

%% To run 1 iteration of the model. 
% See this file for step-by-step model dynamics and learning.
edit HVCIter
%% To get paper figures:
%% Alternating Differentiation
% Code to generate Figure 5 a-f, which shows alternating seed neuron 
% differentiation, from subsong through protosyllable stage through 
% splitting.

edit AlternatingDifferentiation  % run this file

% relies on: 

% HVCIter                   % runs 1 iteration of the model.  See this file for 
                            % step-by-step model dynamics and learning.
                    
% plotSubsong               % plotting function for subsong network diagram and
                            % raster
                    
% plotHVCnet                % plot network diagram

% plotAlternating           % plot raster

% findLatency               % called by plotting functions, tests what neurons
                            % participate in each syllable, and at what latencies
                    
%% Bout onset differentiation
% Code to generate Figure EDF 10 a-d, which shows bout onset differentiation

edit BoutOnsetDifferentiation    % run this file

% relies on: 

% HVCIter                   % runs 1 iteration of the model.  See this file for 
                            % step-by-step model dynamics and learning.
                    
% plotHVCnet_boutOnset      % plotting function for bout onset network diagram 
                            % and raster
                    
% findLatency               % called by plotting functions, tests what neurons
                            % participate in each syllable, and at what latencies

%% latency distribution over development    

edit SigLatDistOverDev.m

%% Bout onset element

edit RunHVC_boutOnsetElement.m

%% Motif learning

edit RunHVC_split_intoThree.m

%% Alternating differentiation movie for supp

edit runHVC_split_movie.m

%% Bout onset movie for supp

edit boutOnsetDifferentiation_movie.m

%% Bout onset movie, sorting weight matrix and multi-dimensional scaling

edit RunHVC_boutOnset_Movies