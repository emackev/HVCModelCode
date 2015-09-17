% Emily Mackevicius 1/14/2015, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

% Code to generate figure EDF10 a-d, which shows bout onset differentiation

clear all; 

%% Bout-onset differentiation: network parameters

% fixed parameters
seed = 1009; 
p.seed = seed;          % seed random number generator
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 500;         % time-steps to simulate -- 
                        %   each time-step is 1 burst duration.
p.pin = .01;            % probability of external stimulation 
                        %   of at least one neuron at any time
k = 10;                 % number of training neurons
p.trainingInd = 1:k;    % index of training neurons
p.beta = .13;           % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .05;            % learning rate parameter
p.epsilon = .14;        % relative strength of heterosynaptic LTD
p.tau = 4;              % time constant of adaptation
gammaStart = .01;       % strength of recurrent inhibition
gammaSplit = .04;       % increased strength of recurrent inhibition 
                        %   to induce splitting
wmaxStart = 1;          % single synapse hard bound
wmaxSplit = 2;          % single synapse hard bound to induce splitting 
                        %  (increased to encourage fewer stronger synapses)
mStart = 5;             % desired number of synapses per neuron 
                        %   (wmax = Wmax/m)
Wmax = mStart*wmaxStart;% soft bound for weights of each neuron
mSplit = Wmax/wmaxSplit;% keep Wmax constant, change m & wmax 
                        %   to induce fewer stronger synapses
HowClamped = 10;        % give training neurons higher threshold
HowOn = 10;             % higher inputs to bout onset training neurons
HowOnPsylStart = HowOn; % inputs to protosyllable training neurons
HowOnPsylSplit = 1;     % decrease input to protosyllable training neurons
                        %   during splitting

% how many iterations to run before plotting
nIterEarly = 5;         % early protosyllable stage
nIterProto = 100;       % end of protosyllable stage
nIterPlotSplit1 = 30;   % number of splitting iterations before plotting 
                        %   intermediate splitting phase
nIterPlotSplit2 = 500;  % total number of splitting iterations

% parameters that change over development
protosyllableStage = [true(1,nIterProto) false(1,nIterPlotSplit2)]; 
splittingStage = [false(1,nIterProto) true(1,nIterPlotSplit2)];
gammas(protosyllableStage) = gammaStart;
gammas(splittingStage) = gammaSplit * sigmf(1:nIterPlotSplit2,[1/200 250]);
wmaxs(protosyllableStage) = wmaxStart;
wmaxs(splittingStage)     = wmaxSplit;
ms(protosyllableStage) = mStart;
ms(splittingStage)     = mSplit;
HowOnPsyl(protosyllableStage) = HowOnPsylStart; 
HowOnPsyl(splittingStage)     = HowOnPsylSplit; 

% params for training inputs
CyclesPerBout = 5; 
bOnOffset = 3; 

%% Bout-onset differentiation: run simulation

% random initial weights
rng(seed);
w = 2*rand(p.n)*Wmax/p.n; 
bOnOffsetVar = [1 randperm(20)]; % variable inter-bout-interval

% learning stages
for t = 1:(nIterProto+nIterPlotSplit2)  
    p.w = w;
    % set parameters that change over development
    p.gamma = gammas(t); 
    p.wmax = wmaxs(t); 
    p.m = ms(t); 
    % Construct input
    Input = -HowClamped*ones(k, p.nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)]; % variable inter-bout-interval
    % initializing
    indPsyl = []; indBstart = []; indOff = []; prevPsylEnd = 1; 
    for i = 1:(p.nsteps/CyclesPerBout/p.trainint)
        istart = (i-1)*CyclesPerBout*p.trainint+1+bOnOffsetVar(i)+bOnOffset; 
        indBstart = [indBstart istart-bOnOffset]; % bout onset times
        indPsyl = [indPsyl ...
            istart istart+p.trainint istart+2*p.trainint]; % 3psyls/bout
        indOff = [indOff ...
            prevPsylEnd:(istart-bOnOffset-1)]; % will clamp all neurons 
                                               % between bouts
        prevPsylEnd = istart+3*p.trainint; % keep track of when bout ends, 
                                           % to clamp neurons between bouts
    end
    indPsyl = indPsyl(indPsyl<=p.nsteps); 
    indBstart = indBstart(indBstart<=p.nsteps);
    Input(1:k/2,indBstart) = HowOn; % input to bout onset neurons
    Input((k/2+1):k,indPsyl) = HowOnPsyl(t); % input to psyl neurons    
    bdyn = double(rand(p.n,p.nsteps)>=(1-p.pin)); % Random activation
    bdyn(1:k,:) = Input; 
    bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
    p.input = bdyn;    
    
    % run one iteration
    [w xdyn] = HVCIter(p);    
    
    % save certain iterations for plotting later
    switch t
        case nIterEarly
            wEarly = w; 
            xdynEarly = xdyn; 
            trainingNeuronsEarly{1}.tind = indBstart+bOnOffset; 
            trainingNeuronsEarly{2}.tind = ...
                setdiff(indPsyl, indBstart+bOnOffset); 
        case nIterProto; 
            wProto = w; 
            xdynProto = xdyn; 
            trainingNeuronsProto{1}.tind = indBstart+bOnOffset; 
            trainingNeuronsProto{2}.tind = ...
                setdiff(indPsyl, indBstart+bOnOffset); 
        case nIterProto + nIterPlotSplit1; 
            wSplit1 = w; 
            xdynSplit1 = xdyn; 
            trainingNeuronsSplit1{1}.tind = indBstart+bOnOffset; 
            trainingNeuronsSplit1{2}.tind = ...
                setdiff(indPsyl, indBstart+bOnOffset); 
        case nIterProto + nIterPlotSplit2;
            wSplit2 = w; 
            xdynSplit2 = xdyn; 
            trainingNeuronsSplit2{1}.tind = indBstart+bOnOffset; 
            trainingNeuronsSplit2{2}.tind = ...
                setdiff(indPsyl, indBstart+bOnOffset); 
    end
end

%% Bout-onset differentiation: plotting parameters

figure(2)

isEPS = 0; 
clf
set(gcf, 'color', ones(1,3));

if isEPS 
    PlottingParams.msize = 8; % change to what is best for EPS figure
    PlottingParams.linewidth = .25;
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultTextFontName', 'Arial')
    PlottingParams.labelFontSize = 7; 
    set(gcf, 'units','centimeters', 'position', [5 5 13.5 6])
else
    PlottingParams.msize = 10;
    PlottingParams.linewidth = .25;
    PlottingParams.labelFontSize = 7; 
end

PlottingParams.SeedColor = [.95 .5 1];
PlottingParams.Syl1Color = [0 0 1]; 
PlottingParams.Syl2Color = [1 0 0]; 
PlottingParams.Syl1BarColor = [0 0 1]; 
PlottingParams.Syl2BarColor = [1 0 0];
PlottingParams.numFontSize = 5; 
PlottingParams.wplotmin = 0; 
PlottingParams.wplotmax = 2; % this should be wmaxSplit
PlottingParams.wprctile = 0; % plot all weights above this percentile. 
PlottingParams.totalPanels = 4; 
PlottingParams.thisPanel = 1; 
PlottingParams.sortby = 'weightMatrix'; 

%% Bout-onset differentiation: plotting early network activity
 
trainingNeuronsEarly{1}.nIDs = 1:k/2;
trainingNeuronsEarly{2}.nIDs = (k/2+1):k;
trainingNeuronsEarly{1}.candLat = (-bOnOffset+1):p.trainint;
trainingNeuronsEarly{2}.candLat =  1:p.trainint; 
trainingNeuronsEarly{1}.thres = 4;
trainingNeuronsEarly{2}.thres = 6;

PlottingParams.thisPanel = 1;
PlottingParams.Hor = 0; 
pp1 = PlottingParams; 
pp1.Syl1BarColor = [1 1 1]; 
pp1.Syl2BarColor = [.5 .5 .5];
plotHVCnet_boutOnset(wEarly, xdynEarly, trainingNeuronsEarly, pp1)
PlottingParams.Hor = 1;

%% Bout-onset differentiation: plotting protosyllable
trainingNeuronsProto{1}.nIDs = 1:k/2;
trainingNeuronsProto{2}.nIDs = (k/2+1):k;
trainingNeuronsProto{1}.candLat = (-bOnOffset+1):p.trainint;
trainingNeuronsProto{2}.candLat =  1:p.trainint; 
trainingNeuronsProto{1}.thres = 4;
trainingNeuronsProto{2}.thres = 6;

PlottingParams.thisPanel = 2;
plotHVCnet_boutOnset(wProto, xdynProto, ...
    trainingNeuronsProto, PlottingParams)

%% Bout-onset differentiation: plotting splitting stages
trainingNeuronsSplit1{1}.nIDs = 1:k/2;
trainingNeuronsSplit1{2}.nIDs = (k/2+1):k;
trainingNeuronsSplit1{1}.candLat = (-bOnOffset+1):p.trainint;
trainingNeuronsSplit1{2}.candLat =  1:p.trainint; 
trainingNeuronsSplit1{1}.thres = 4;
trainingNeuronsSplit1{2}.thres = 6;

PlottingParams.thisPanel = 3;
plotHVCnet_boutOnset(wSplit1, xdynSplit1, ...
    trainingNeuronsSplit1, PlottingParams)

trainingNeuronsSplit2{1}.nIDs = 1:k/2;
trainingNeuronsSplit2{2}.nIDs = (k/2+1):k;
trainingNeuronsSplit2{1}.candLat = (-bOnOffset+1):p.trainint;
trainingNeuronsSplit2{2}.candLat =  1:p.trainint; 
trainingNeuronsSplit2{1}.thres = 4;
trainingNeuronsSplit2{2}.thres = 6;

PlottingParams.thisPanel = 4;
plotHVCnet_boutOnset(wSplit2, xdynSplit2, ...
    trainingNeuronsSplit2, PlottingParams)

%% Bout-onset differentiation: exporting
if isEPS
    cd('Z:\Fee_lab\Papers\HVC_differentiation\Figures\EPS_files');
    export_fig(2,'Figure5h.eps','-transparent','-eps','-painters');
else
    figw = 6;
    figh = 4;
    set(gcf, 'color', [1 1 1],...
        'papersize', [figw figh], 'paperposition', [0 0 figw*.9 figh])
    %print -dmeta -r150
end
