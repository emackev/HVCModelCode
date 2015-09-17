function HVCModelCode
% Emily Mackevicius 7/18/2015, based on Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 
%% *Alternating differentiation*
% From subsong through protosyllable stage through splitting, 
% to generate figure 5 a-d
% Calls HVCIter to step through one iteration of the model

%% Alternating differentiation: network parameters

% fixed parameters
seed = 9038;
p.seed = seed;          % seed random number generator
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 100;         % time-steps to simulate -- each time-step is 1 
                        %   burst duration.
nstepsSubsong = 1000;   % time-steps to simulate for subsong stage
p.pin = .01;            % probability of external stimulation of at least 
                        %   one neuron at any time
k = 10;                 % number of training neurons
p.trainingInd = 1:k;    % index of training neurons
p.beta = .115;          % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .025;           % learning rate parameter
p.epsilon = .2;         % relative strength of heterosynaptic LTD
p.tau = 4;              % time constant of adaptation
gammaStart= .01;        % strength of recurrent inhibition
gammaSplit =.18;        % increased strength of recurrent inhibition to 
                        %   induce splitting
wmaxStart = 1;          % single synapse hard bound
wmaxSplit = 2;          % single synapse hard bound to induce splitting 
                        %  (increased to encourage fewer stronger synapses)
mStart = 10;            % desired number of synapses per neuron 
                        %   (wmax = Wmax/m)
Wmax = mStart*wmaxStart;% soft bound for weights of each neuron
mSplit = Wmax/wmaxSplit;% keep Wmax constant, change m & wmax to induce 
                        %   fewer stronger synapses
HowClamped = 10;        % give training neurons higher threshold
HowOn = 10;             % higher inputs to training neurons

% how many iterations to run before plotting
nIterProto = 500;       % end of protosyllable stage
nIterPlotSplit1 = 492;  % number of splitting iterations before plotting 
                        %   intermediate splitting phase
nIterPlotSplit2 = 2000; % total number of splitting iterations

% parameters that change over development
protosyllableStage = [true(1,nIterProto) false(1,nIterPlotSplit2)]; 
splittingStage = [false(1,nIterProto) true(1,nIterPlotSplit2)];
gammas(protosyllableStage) = gammaStart;
gammas(splittingStage) = gammaSplit*sigmf(1:nIterPlotSplit2,[1/200 500]);
wmaxs(protosyllableStage) = wmaxStart;
wmaxs(splittingStage)     = wmaxSplit;
ms(protosyllableStage) = mStart;
ms(splittingStage)     = mSplit;

%Subsong Inputs
rng(seed)
isOnset = rand(1,nstepsSubsong)>.9; 
Input =-HowClamped*ones(k, nstepsSubsong); % clamp training neurons 
                        % (effectively giving them higher threshold)
Input(:,isOnset) = HowOn; 
bdyn = double(rand(p.n,nstepsSubsong)>=(1-p.pin)); % Random activation
bdyn(1:k,:) = Input; 
subsongInput = bdyn; 

%Protosyllable inputs
PsylInput = -HowClamped*ones(k, p.nsteps); %clamp training neurons 
                        % (effectively giving them higher threshold)
PsylInput(:,mod(1:p.nsteps,p.trainint)==1) = HowOn; % rhythmic activation 
                                                    %   of training neurons

%Alternating Inputs
AltInput =-HowClamped*ones(k, p.nsteps); 
AltInput(1:k/2,mod(1:p.nsteps,2*p.trainint)==1) = HowOn; 
AltInput((k/2+1):k,mod(1:p.nsteps,2*p.trainint)==p.trainint+1) = HowOn; 
                    % alternating rhythmic activation of training neurons

%% Alternating differentiation: run simulation

% random initial weights
rng(seed);
w0 = 2*rand(p.n)*Wmax/p.n;

% subsong stage
pSubsong = p; 
pSubsong.gamma = gammas(1); 
pSubsong.wmax = wmaxs(1); 
pSubsong.m = ms(1); 
pSubsong.eta = 0; 
pSubsong.epsilon = 0; 
pSubsong.nsteps = nstepsSubsong; 
pSubsong.w = w0; 
pSubsong.input = subsongInput;

% Run subsong network
[wSubsong xdynSubsong] = HVCIter(pSubsong);
w = wSubsong;

% learning stages
for t = 1:(nIterProto+nIterPlotSplit2)  
    p.w = w;
    % set parameters that change over development
    p.gamma = gammas(t); 
    p.wmax = wmaxs(t); 
    p.m = ms(t); 
    % Construct input
    bdyn = double(rand(p.n,p.nsteps)>=(1-p.pin)); % Random activation
    bdyn(1:k,:) = protosyllableStage(t)*PsylInput+ ...
        splittingStage(t)*AltInput; % drive to seed neurons
    p.input = bdyn;
    % run one iteration
    [w xdyn] = HVCIter(p);    
    % save certain iterations for plotting later
    switch t
        case nIterProto; 
            wProto = w; 
            xdynProto = xdyn; 
        case nIterProto + nIterPlotSplit1; 
            wSplit1 = w; 
            xdynSplit1 = xdyn; 
        case nIterProto + nIterPlotSplit2;
            wSplit2 = w; 
            xdynSplit2 = xdyn; 
    end
end

%% Alternating differentiation: plotting parameters
figure(1)
isEPS = 0; 
clf
set(gcf, 'color', ones(1,3));

if isEPS 
    PlottingParams.msize = 8; 
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
PlottingParams.Syl1Color = [1 0 0]; 
PlottingParams.Syl2Color = [0 0 1]; 
PlottingParams.ProtoSylColor = [0 0 0]; 
PlottingParams.ProtoSylBarColor = [.5 .5 .5];
PlottingParams.SubsongSylColor = [0 0 0]; 
PlottingParams.SubsongBarColor = [1 1 1]; 
PlottingParams.numFontSize = 5; 
PlottingParams.wplotmin = 0; 
PlottingParams.wplotmax = 2; % this should be wmaxSplit
PlottingParams.wprctile = 0; % plot all weights above this percentile.  
                             %  If nonzero, ignores wplotmin, wplotmax
PlottingParams.wperneuron = 6; % max outgoing weights plotted
PlottingParams.wperneuronIn = 9; % min incoming weights plotted
PlottingParams.totalPanels = 4; 
nplots = 4;
bottom = .1; 
height = .45; 
scale = .005; 
spacing = .75/(2*nplots); 

%% Alternating differentiation: plotting subsong
trainingNeuronsSubsong{1}.nIDs = 1:k;
trainingNeuronsSubsong{1}.tind = find(isOnset);
trainingNeuronsSubsong{1}.candLat = 1:2*p.trainint; 
trainingNeuronsSubsong{1}.thres = 12; % criteria for participation during 
                                      %  subsong (thres from testLatSig -- 
                                      %  must fire at consistent latency 
                                      %  more than 12 times in the bout of 
                                      %  ~100 syllables to count as 
                                      %  participating)

PlottingParams.thisPanel = 1; 
PlottingParams.Hor = 1;
plotSubsong(wSubsong, xdynSubsong, trainingNeuronsSubsong, PlottingParams)

%% Alternating differentiation: plotting protosylable
trainingNeuronsPsyl{1}.nIDs = 1:k;
trainingNeuronsPsyl{2}.nIDs = 1:k;
trainingNeuronsPsyl{1}.tind = find(mod(1:p.nsteps, p.trainint)==1);
trainingNeuronsPsyl{2}.tind = find(mod(1:p.nsteps, p.trainint)==1);
trainingNeuronsPsyl{1}.candLat = 1:p.trainint; 
trainingNeuronsPsyl{2}.candLat = 1:p.trainint; 
trainingNeuronsPsyl{1}.thres = 4; % criteria for participation during 
                                  %  protosyllable stage (thres from 
                                  %  testLatSig -- must fire at consistent 
                                  %  latency more than 4 times in the bout 
                                  %  of 10 syllables to count as 
                                  %  participating)
trainingNeuronsPsyl{2}.thres = 4; 

PlottingParams.thisPanel = 2; 
subplot('position', ...
    [PlottingParams.thisPanel/nplots-.9/nplots, .6, .9/nplots, .35])
plotHVCnet(wProto,xdynProto,p.trainint,trainingNeuronsPsyl,PlottingParams)
set(gca, 'color', 'none');
PlottingParams.axesPosition = ...
    [PlottingParams.thisPanel/nplots-2*spacing, bottom, 40*scale, height];
plotAlternating(wProto, xdynProto, p.trainint, ...
    trainingNeuronsPsyl, PlottingParams)
set(gca, 'color', 'none')

%% Alternating differentiation: plotting splitting stages
trainingNeuronsAlt{1}.nIDs = 1:k/2;
trainingNeuronsAlt{2}.nIDs = (k/2+1):k;
trainingNeuronsAlt{1}.tind = find(mod(1:p.nsteps, 2*p.trainint)==1);
trainingNeuronsAlt{2}.tind = find(mod(1:p.nsteps, ...
                                2*p.trainint)==p.trainint+1);
trainingNeuronsAlt{1}.candLat = 1:p.trainint; 
trainingNeuronsAlt{2}.candLat = 1:p.trainint; 
trainingNeuronsAlt{1}.thres = 2; % criteria for participation during 
                                 %  splitting stage (thres from testLatSig 
                                 %  -- must fire at consistent latency 
                                 %  more than 2 times in the bout of 5 
                                 %  syllables (of each type) to count 
                                 %  as participating)
trainingNeuronsAlt{2}.thres = 2; 

PlottingParams.thisPanel = 3;
subplot('position', ...
    [PlottingParams.thisPanel/nplots-.9/nplots, .6, .9/nplots, .35])
plotHVCnet(wSplit1,xdynSplit1,p.trainint,trainingNeuronsAlt,PlottingParams)
set(gca, 'color', 'none');
PlottingParams.axesPosition = ...
    [PlottingParams.thisPanel/nplots-2*spacing, bottom, 40*scale, height];
plotAlternating(wSplit1, xdynSplit1, p.trainint, ...
    trainingNeuronsAlt, PlottingParams)
set(gca, 'color', 'none')

PlottingParams.thisPanel = 4;
subplot('position', ...
    [PlottingParams.thisPanel/nplots-.9/nplots, .6, .9/nplots, .35])
plotHVCnet(wSplit2,xdynSplit2,p.trainint,trainingNeuronsAlt,PlottingParams)
set(gca, 'color', 'none');
PlottingParams.axesPosition = ...
    [PlottingParams.thisPanel/nplots-2*spacing, bottom, 40*scale, height];
plotAlternating(wSplit2, xdynSplit2, p.trainint, ...
    trainingNeuronsAlt, PlottingParams)
set(gca, 'color', 'none')

%% Alternating differentiation: exporting
if isEPS
    cd('Z:\Fee_lab\Papers\HVC_differentiation\Figures\EPS_files');
    export_fig(1,'Figure5a.eps','-transparent','-eps','-painters');
else
    %figure parameters, exporting
    figw = 6;
    figh = 2;
    set(gcf, 'color', [1 1 1],...
        'papersize', [figw figh], 'paperposition', [0 0 figw*.9 figh])
    % print -dmeta -r150
end

%% *Bout-onset differentiation* 

% Code to generate figure 5 j-m, which shows bout onset differentiation

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

end

%% *HVCIter function*

function [w xdyn] = HVCIter(p)
% Runs one iteration of the simulation.  p is a structure of parameters.

% redefining params that are used often outside the for loop
nsteps = p.nsteps;
n = p.n;
m = p.m;
w = p.w; 
wmax = p.wmax;
Wmax = wmax*m; 
eta = p.eta;
bdyn=p.input;

% initializing variables
xdyn=zeros(n,nsteps);
oldx = zeros(n,1);
oldy = zeros(n,1);

for t = 1:nsteps
    % Adaptation
    y = oldy + 1/p.tau*(-oldy+oldx);
    Aadapt = p.alpha * y; % adaptation
    
    % Net feedforward input.  
    B = bdyn(:,t); % external input
    AE = w*oldx; % excitatory input 
    AIff = p.beta * sum(oldx); % feed forward inhibition
    Anet = AE - AIff - Aadapt + B; % net feed forward input
    Anet(Anet < 0) = 0; % rectify
    
    % recurrent inhibition
    AIrec = p.gamma * sum(Anet); 
    
    % binary output
    x = (Anet - AIrec) > 0; 
    
    % STDP rule (Fiete et al 2010)
    dw_STDP = eta.*(x*(oldx)'-(oldx)*x');
    
    % Hetersynaptic penalty (Fiete et al 2010)
    dw_hLTDpre = ...
        eta*ones(n,1)*max(0, sum(w+dw_STDP,1)-Wmax);  % Weights leaving 
                                                      %  cells (pre)
    dw_hLTDpost = ...
        eta*max(0, sum(w+dw_STDP,2)-Wmax)*ones(1,n);  % Weights onto 
                                                      %  cells (post)
    
    % Update weights
    if eta>0
        dwtotal = dw_STDP-p.epsilon*(dw_hLTDpre+dw_hLTDpost);
        w = w + dwtotal;
        w(w > wmax) = wmax; % hard limit on strength of a single synapse
        w(w < 0) = 0; % weights cannot be negative
        w = w.*(~eye(p.n)); % clamp diagonal
    end
    
    oldx = double(x);
    oldy = y;
    xdyn(:,t)=x;
end
end

%% Plotting functions: *findLatency*
function Latency = findLatency(xsort, trainingNeurons)
% Calculates the mode latency of each neuron
% w: weight matrix
% xdyn: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing 
%   neuron and time indices for each syllable type
% PlottingParams: sets linewidth, etc.  

nsteps = size(xsort,2); 
n = size(xsort,1); 

% iterate over candidate latencies
for syli = 1:length(trainingNeurons)
    nFired = zeros(n,length(trainingNeurons{syli}.candLat));
    for lati = 1:length(trainingNeurons{syli}.candLat); 
        tInds = zeros(1,nsteps); 
        tInds(min(nsteps, ...
            trainingNeurons{syli}.tind + ...
            trainingNeurons{syli}.candLat(lati))-1) = 1; 
        nFired(:,lati) = ...
            sum(bsxfun(@times, xsort, tInds),2); % number of times each 
                                                 % neuron fired at this 
                                                 % latency
    end
    for ni = 1:n
        [~, Latency{syli}.mode(ni)] = max(nFired(ni,:)); 
        Latency{syli}.mode(ni) = ...
            trainingNeurons{syli}.candLat(Latency{syli}.mode(ni)); 
        Latency{syli}.FireDur(ni) = ...
            max(nFired(ni,:))>trainingNeurons{syli}.thres;  
    end
end
end


%% Plotting functions: *plotAlternating* 
function plotAlternating(w, xsort, m, trainingNeurons, PlottingParams)
% Makes network activity plot, called by AlternatingDifferentiation
% w: weight matrix
% xsort: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing 
%   neuron and time indices for each syllable type
% PlottingParams: sets linewidth, etc.  

Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
ProtoSylColor = PlottingParams.ProtoSylColor;
numFontSize = PlottingParams.numFontSize;
labelFontSize = PlottingParams.labelFontSize;

Latency = findLatency(xsort, trainingNeurons);

% plotting the mode latency for each syll type
xplot = zeros(size(w,1),2*m);
for ni = 1:size(w,1) 
    if Latency{1}.FireDur(ni) & ~isnan(Latency{1}.mode(ni))
        xplot(ni,Latency{1}.mode(ni)) = 1; 
    end
    if Latency{2}.FireDur(ni) & ~isnan(Latency{2}.mode(ni))
        xplot(ni,Latency{2}.mode(ni)+m) = 1; 
    end
end

% keep track of which neurons participated in each syllable
FireDur1 = Latency{1}.FireDur; 
FireDur2= Latency{2}.FireDur; 

% classify neurons as specific or shared
Specific1 = FireDur1&~FireDur2;
Specific2 = FireDur2&~FireDur1; 
Shared = (FireDur1&FireDur2); 
indshared = find(Shared);

Red = trainingNeurons{1}.nIDs;
Green = trainingNeurons{2}.nIDs; 

if issame(Red,Green)
    IsProto = zeros(1,length(xplot)); IsProto(Red) = 1;
    Red = [];
    Green = []; 
else
    IsProto = zeros(1,length(xplot));
end
IsTrain1 = zeros(1,length(xplot)); IsTrain1(Red) = 1; 
IsTrain2 = zeros(1,length(xplot)); IsTrain2(Green) = 1;


if length(Green)==0
    tmp= xplot>0;
    tmpXplot = tmp(:,1:(size(xplot,2)/2))+tmp(:,(size(xplot,2)/2+1):end);
else
    tmpXplot = xplot>0;
end
[~,sortind] = sortrows(tmpXplot); 
xplot = xplot(flipud(sortind),:); % from early (top) to late (bottom)
w = w(flipud(sortind),flipud(sortind));
IsTrain1 = IsTrain1(flipud(sortind)); 
IsTrain2 = IsTrain2(flipud(sortind)); 
Specific1 = Specific1(flipud(sortind)); 
Specific2 = Specific2(flipud(sortind)); 
IsProto = IsProto(flipud(sortind)); 

% if differentiated, sort shared neurons first, then specific neurons
if length(Green)>0
    sharedind = (sum(xplot,2)>=2)&(sum(xplot,2)<8);
else
    sharedind = zeros(1,size(xplot,1)); 
end
rest = find(~sharedind); 
xplotall = xplot([find(sharedind); (rest)],:);
IsTrain1 = IsTrain1([find(sharedind); (rest)]); 
IsTrain2 = IsTrain2([find(sharedind); (rest)]); 
Specific1 = Specific1([find(sharedind); (rest)]); 
Specific2 = Specific2([find(sharedind); (rest)]); 
IsProto = IsProto([find(sharedind); (rest)]); 

% plotting the activity for the two syll types
for i = 1:2
    axesPos = PlottingParams.axesPosition;
    axesPos(1) = axesPos(1)+(i-1)*axesPos(3)/2; 
    axesPos(3) = axesPos(3)*.3;
    subplot('position', axesPos); 
    xplot = ...
        xplotall(:,(1:(size(xplotall,2)/2))+(i-1)*(size(xplotall,2)/2));
    tplot = (1:(size(xplot,2)))*10; % assuming each bin is 10ms

    tOffset = 0; 
    for j=1:size(xplot,2) % for all the time steps
        Idx = ...
            find(xplot(1:end-1,j)>0); % find the indices of active neurons    
        if ~isempty(Idx)
            for k=1:length(Idx) % for all the active neurons
                Color = IsProto(Idx(k))*PlottingParams.ProtoSylColor + ...
                    Specific1(Idx(k))*PlottingParams.Syl1Color + ...
                    Specific2(Idx(k))*PlottingParams.Syl2Color;
                h = patch(10*([j-1,j,j,j-1]+tOffset),...
                    [Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],...
                    Color,'edgecolor','none');
            end  
        end
    end
    hold on
    
    % plot line between each syl type
    if length(Green)>0
        train2 = find(IsTrain2); train2 = train2(1)-1; 
        plot([0 size(xplot,2)*10], [train2 train2], 'k', ...
            'linewidth', PlottingParams.linewidth) 
    end
    
    % plot line between shared and specific neurons
    if length(rest)>0 & sum(sharedind)>0 & length(Green)>0
        plot([0 size(xplot,2)*10], [sum(sharedind) sum(sharedind)],...
            'k', 'linewidth', PlottingParams.linewidth) 
    end
    
    % plot colored bars above each syllable
    if length(Green) == 0
        patch([0 90 90 0],[-4 -4 -2 -2],PlottingParams.ProtoSylBarColor); 
        text(40,-10,'\alpha','fontsize',7); 
    elseif i == 1
        patch([0 90 90 0],[-4 -4 -2 -2],Syl1Color); 
        text(40,-10,'\beta','fontsize',7) 
    else
        patch([0 90 90 0],[-4 -4 -2 -2],Syl2Color); 
        text(40,-10,'\gamma','fontsize',7); 
    end
    
    % plotting parameters
    ylim([-5 size(xplot,1)-1]);  
    box off
    set(gca, 'ydir', 'reverse','tickdir','out',...
        'ticklength',[0.015 0.015], 'color', 'none', ...
        'xtick', 0:50:100,'fontsize', numFontSize,'tickdir','out');
    xlim([-2 100]);
    if PlottingParams.thisPanel==1
        ylabel('Neuron #', 'fontsize', labelFontSize,'fontname','arial');
        set(gca,'ytick',0:20:100,'fontsize',numFontSize)
    else
        set(gca,'ytick',0:20:100,'yticklabel', {})
    end
end
end

%% Plotting functions: *plotHVCnet*
function plotHVCnet(w, xdyn, m, trainingNeurons, PlottingParams)
% Makes network diagram for alternating differentiation
% w: weight matrix
% xdyn: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing 
%   neuron and time indices for each training neuron type
% PlottingParams: sets linewidth, etc. 

msize = PlottingParams.msize;
linewidth = PlottingParams.linewidth;
Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
ProtoSylColor = PlottingParams.ProtoSylColor;

Latency = findLatency(xdyn, trainingNeurons);

% first exclude all neurons that don't fire at a consistent phase
cla; hold on
x = zeros(1,size(w,1));
y = zeros(1,size(w,1));
for ni = 1:size(w,1)
    % if it fired during either syll
    if Latency{1}.FireDur(ni)|Latency{2}.FireDur(ni) 
        % if it fired during both sylls    
        if (Latency{1}.FireDur(ni)&Latency{2}.FireDur(ni)) 
            % if fired during both sylls at same phase    
            if (Latency{1}.mode(ni)==Latency{2}.mode(ni)) 
                x(ni) = Latency{1}.mode(ni);
            else % exclude from plot if different phases for both sylls
                x(ni) = NaN;
            end
        elseif Latency{1}.FireDur(ni) % if it fired during syll 1 
            x(ni) = Latency{1}.mode(ni);
        else % fired during syll 2 only
            x(ni) = Latency{2}.mode(ni);
        end
    else % if it fired during neither syll
        x(ni) = NaN;
    end
end
indkeep = ~isnan(x); 
y = y(indkeep); 
w = w(indkeep,indkeep); 
xdyn = xdyn(indkeep,:); 
x = x(indkeep); 
ux = unique(x); 

% keep track of training neuron and syl time indices
trainingset1 = trainingNeurons{1}.nIDs;
trainingset2 = trainingNeurons{2}.nIDs; 
x(trainingset1) = 1; 
x(trainingset2) = 1; 
tind1 = (trainingNeurons{1}.tind);
tind2 = (trainingNeurons{2}.tind);

% keep track of which neurons participated in each syllable
FireDur1 = Latency{1}.FireDur(indkeep); 
FireDur2= Latency{2}.FireDur(indkeep); 

% classify neurons as specific or shared
Specific1 = FireDur1&~FireDur2;
Specific2 = FireDur2&~FireDur1; 
Shared = (FireDur1&FireDur2); 
indshared = find(Shared);

% calculate the incoming weights from specific neurons of each type, to
%  determine sorting in y axis and color
c1 = zeros(1,length(x));
c2 = zeros(1,length(x));
for ni = 1:size(w,1)
    tmp = find(xdyn(ni,:)); 
    if sum(w(ni,:))>0
        c1(ni) = sum(w(ni,Specific1))/sum(w(ni,:));
        c2(ni) = sum(w(ni,Specific2))/sum(w(ni,:));
    end
    y(ni) = c1(ni)-c2(ni);
end

% for each latency (x), sort along y, with small gap between shared and
%  specific neurons
y1 = zeros(1,size(w,1)); 
for ui = 1:length(ux)
    indshared = (x==ux(ui))&Shared;
    ind1 = (x==ux(ui))&Specific1; 
    ind2 = (x==ux(ui))&Specific2;
    [~,y1(indshared)] = sort(y(indshared));
    tocentershared = 1+(numel(find(indshared))-1)/2;
    y1(indshared) = y1(indshared)-tocentershared;
    [~,y1(ind1)] = sort(y(ind1));
    y1(ind1) = y1(ind1) + (numel(find(indshared)))/2;
    [~,y1(ind2)] = sort(y(ind2));
    y1(ind2) = y1(ind2) -numel(find(ind2))-1- (numel(find(indshared)))/2;
end

cla; hold on

% keep only feedforward part of weight matrix
wplot = w; 
n = size(wplot,1); 
for i = 1:n
    for j = 1:n
        ff = x(i)<x(j);
        longrange = abs(x(i)-x(j))>2; 
        if (~ff) | longrange
            wplot(j,i) = 0; 
        end
    end
end

% Color weights white to black between wplotmin and wplotmax
wplot = wplot-PlottingParams.wplotmin; 
wplot(wplot<0) = 0;
wplot = wplot/(PlottingParams.wplotmax-PlottingParams.wplotmin);
wplot(wplot<prctile(wplot(:), PlottingParams.wprctile)) = 0; 
wplotold = wplot; 
n = size(wplot,1); 
for i = 1:n
    [~,ind] = sort(wplot(:,i), 'descend'); 
    indplot = zeros(n,1); 
    indplot(ind(1:min(PlottingParams.wperneuron,length(ind)))) = 1;   
    wplot(~indplot,i) = 0; 
end
for i = 1:n 
    if sum(wplot(i,:)>0)<PlottingParams.wperneuronIn
        [~,ind] = sort(wplotold(i,:), 'descend'); 
        wplot(i,ind(1:min(PlottingParams.wperneuron,length(ind)))) = ...
            wplotold(i,ind(1:min(PlottingParams.wperneuron,length(ind))));
    end
end

% jitter a little in x and y, so it doesn't look like a grid, but
%  don't jitter seed neurons
jitter = .1; 
Seed0 = randn(1,300);
indJitter = setdiff(1:length(x), union(trainingset1, trainingset2)); 
x(indJitter)= x(indJitter)+jitter*Seed0(1:length(x(indJitter)));
y1(indJitter) = ...
    y1(indJitter)+...
    jitter*Seed0((length(x(indJitter))+1):(2*length(x(indJitter))));

% plot w in order from weakest to strongest, so darker lines are on top
js = repmat((1:n)',1,n); 
is = repmat((1:n),n,1); 
isVec = is(:);
jsVec = js(:); 
wVec = wplot(:); 
[wSort,indSort] = sort(wVec, 'ascend'); 
nplotted = zeros(1,n); 
for k = 1:length(wSort)
    i = isVec(indSort(k)); 
    j = jsVec(indSort(k)); 
    if wplot(j,i)>0
        ff = x(i)<=x(j); 
        longrange = abs(x(i)-x(j))>2; 
        loopback = (round(x(i))==round(max(x)))&...
            (round(x(j))==round(min(x)));
        if (ff & ~longrange)%|loopback
            C = ones(1,3)-wplot(j,i)*ones(1,3);
            plot([x(i), x(j)], [y1(i),y1(j)], ...
                'color', C, 'linewidth', linewidth)
        end
    end
end

% color each neuron based on its relative input from each syllable type
for pli = 1:length(x)
    tmpC = c1(pli)'/(max(c1)+eps)*Syl1Color+c2(pli)'/(max(c2)+eps)*Syl2Color; 
    tmpC = tmpC/(max(tmpC)+eps); % normalize so colors are bright
    if Shared(pli)
        tmpC = zeros(1,3); 
    end
    if Specific1(pli)
        tmpC = Syl1Color; 
    end
    if Specific2(pli)
        tmpC = Syl2Color; 
    end
    plot(x(pli),y1(pli), 'marker', '.', 'color', tmpC, 'markersize', msize)
end

% plot training neurons in given colors
if sum(Specific1)>0
    plot(x(trainingset1),y1(trainingset1), ...
        '.', 'markersize', msize, 'color', Syl1Color)
    plot(x(trainingset2),y1(trainingset2), ...
        '.', 'markersize', msize, 'color', Syl2Color)
else
    plot(x([trainingset1 trainingset2]),y1([trainingset1 trainingset2]),...
        '.', 'markersize', msize, 'color', ProtoSylColor)
end

% plot rectangle for syl1 seed neurons
rx = 1-.5;
ry = min(y1(trainingset1))-.5;
rw = 1;
rh = max(y1(trainingset1)) - min(y1(trainingset1))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor', 'none',...
    'LineStyle', '-', 'LineWidth', .5, ...
    'EdgeColor', PlottingParams.SeedColor, ...
    'curvature', [.98 .1])

% plot rectangle for syl2 seed neurons
rx = 1-.5;
ry = min(y1(trainingset2))-.5;
rw = 1;
rh = max(y1(trainingset2)) - min(y1(trainingset2))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor', 'none',...
    'LineStyle', '-', 'LineWidth', .5, ...
    'EdgeColor', PlottingParams.SeedColor, ...
    'curvature', [.98 .1])

axis tight; axis off; 
xlim([-.5 m+.5]);
ylim([min(y1)-1 max(y1)+1])
set(gca, 'color', 'none')
end

%% Plotting functions: *plotHVCnet_boutOnset*
function plotHVCnet_boutOnset(w, xdyn, trainingNeurons, PlottingParams)
% Makes network diagram and raster plots, called by RunHVC_boutOnset_net
% w: weight matrix
% xdyn: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing 
%   neuron and time indices for each syllable type
% PlottingParams: sets linewidth, etc.  See RunHVC_split

%plotting parameters%
msize = PlottingParams.msize;
linewidth = PlottingParams.linewidth;
Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
numFontSize = PlottingParams.numFontSize;
labelFontSize = PlottingParams.labelFontSize;
nplots = PlottingParams.totalPanels;
ploti = PlottingParams.thisPanel;

%Network diagram%
subplot('position', [ploti/nplots-.9/nplots, .56, .9/nplots, .44])
cla; hold on

% calculate latency of each neuron
Latency = findLatency(xdyn, trainingNeurons);

% first double plot all neurons that don't fire at a consistent phase
nsteps = size(xdyn,2);
n = size(xdyn,1);
ntot = n;
x = zeros(1,n);
y = zeros(1,n);
trainingset1 = trainingNeurons{1}.nIDs;
trainingset2 = trainingNeurons{2}.nIDs;
indDoubled = [];
cDoub = n+1;
for ni = 1:n
    if length(intersect(trainingset1,ni))>0 % if it's a training neuron
        x(ni) = trainingNeurons{1}.candLat(1);
    elseif length(intersect(trainingset2,ni))>0 % if it's a training neuron
        x(ni) = trainingNeurons{2}.candLat(1);
    else
        % if it fired during either syll
        if Latency{1}.FireDur(ni)|Latency{2}.FireDur(ni) 
            % if it fired during both sylls
            if (Latency{1}.FireDur(ni)&Latency{2}.FireDur(ni)) 
                % if fired during both sylls at same phase
                if (Latency{1}.mode(ni)==Latency{2}.mode(ni)) 
                    x(ni) = Latency{1}.mode(ni);
                else % double plot if different phases for both sylls
                    x(ni) = Latency{1}.mode(ni);
                    x(cDoub) = Latency{2}.mode(ni);
                    indDoubled = [indDoubled ni];
                    cDoub = cDoub+1;
                end
            elseif Latency{1}.FireDur(ni) % if it fired during syll 1 only
                x(ni) = Latency{1}.mode(ni);
            else % fired during syll 2 only
                x(ni) = Latency{2}.mode(ni);
            end
        else % if it fired during neither syll
            x(ni) = NaN;
        end
    end
end
indkeep = [find(~isnan(x(1:n))) indDoubled];
y = y(indkeep);
w = w(indkeep,indkeep);
xdyn = xdyn(indkeep,:);
x = x(indkeep);
ux = unique(x);

% keep track of which neurons participated in each syllable
FireDur1 = Latency{1}.FireDur(indkeep);
FireDur2= Latency{2}.FireDur(indkeep);

% classify neurons as specific or shared
Specific1 = FireDur1&~FireDur2;
Specific2 = FireDur2&~FireDur1;
Shared = (FireDur1&FireDur2);
indshared = find(Shared);

% calculate the incoming weights from specific neurons of each type, to
%  determine sorting in y axis and color
c1 = zeros(1,length(x));
c2 = zeros(1,length(x));
for ni = 1:size(w,1)
    tmp = find(xdyn(ni,:));
    if sum(w(ni,:))>0
        c1(ni) = sum(w(ni,Specific1))/sum(w(ni,:));
        c2(ni) = sum(w(ni,Specific2))/sum(w(ni,:));
    end
    y(ni) = c1(ni)-c2(ni);
end

% for each latency (x), sort along y, with small gap between shared and
%  specific neurons
y1 = zeros(1,size(w,1));
for ui = 1:length(ux)
    indshared = (x==ux(ui))&Shared;
    ind1 = (x==ux(ui))&Specific1;
    ind2 = (x==ux(ui))&Specific2;
    [~,y1(indshared)] = sort(y(indshared));
    tocentershared = 1+(numel(find(indshared))-1)/2;
    y1(indshared) = y1(indshared)-tocentershared;
    [~,y1(ind1)] = sort(y(ind1));
    y1(ind1) = y1(ind1) + (numel(find(indshared)))/2;
    [~,y1(ind2)] = sort(y(ind2));
    y1(ind2) = y1(ind2) -numel(find(ind2))-1-(numel(find(indshared)))/2;
end

% Color weights white to black between wplotmin and wplotmax
wplot = w;
wplot = w-PlottingParams.wplotmin;
wplot(wplot<0) = 0;
wplot = wplot/(PlottingParams.wplotmax-PlottingParams.wplotmin);
wplot(wplot<prctile(wplot(:), PlottingParams.wprctile)) = 0;

% jitter a little in x and y, so it doesn't look like a grid, 
%  but don't jitter seed neurons
jitter = .1;
Seed0 = randn(1,500);
indJitter = setdiff(1:length(x), union(trainingset1, trainingset2)); 
x(indJitter)= x(indJitter)+jitter*Seed0(1:length(x(indJitter)));
y1(indJitter) = y1(indJitter)+...
    jitter*Seed0((length(x(indJitter))+1):(2*length(x(indJitter))));

% plot w in order from weakest to strongest, so darker lines are on top
n = size(wplot,1);
js = repmat((1:n)',1,n);
is = repmat((1:n),n,1);
isVec = is(:);
jsVec = js(:);
wVec = wplot(:);
[wSort,indSort] = sort(wVec, 'ascend');
for k = 1:length(wSort)
    i = isVec(indSort(k));
    j = jsVec(indSort(k));
    if wplot(j,i)>0
        ff = x(i)<=x(j);
        longrange = abs(x(i)-x(j))>2;
        loopback = (round(x(i))==round(max(x)))&...
            (round(x(j))==round(min(x)));
        if ff & ~longrange%|loopback
            C = ones(1,3)-wplot(j,i)*ones(1,3);
            plot([x(i), x(j)], [y1(i),y1(j)], ...
                'color', C, 'linewidth', linewidth)
        end
    end
end

% color each neuron based on its relative input from each syllable type
for pli = 1:length(x)
    if Shared(pli)
        tmpC = zeros(1,3);
    end
    if Specific1(pli)
        tmpC = Syl1Color;
    end
    if Specific2(pli)
        tmpC = Syl2Color;
    end
    plot(x(pli),y1(pli), 'marker', '.', 'color', tmpC, 'markersize', msize)
end

% plot rectangle for syl1 seed neurons
rx = mean(x(trainingset1))-.5;
ry = min(y1(trainingset1))-.5;
rw = 1;
rh = max(y1(trainingset1)) - min(y1(trainingset1))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor',  'none', ...
    'LineStyle', '-', 'LineWidth', .5, ...
    'EdgeColor', PlottingParams.SeedColor,...
    'curvature', [.98 .1])

% plot rectangle for syl2 seed neurons
rx = mean(x(trainingset2))-.5;
ry = min(y1(trainingset2))-.5;
rw = 1;
rh = max(y1(trainingset2)) - min(y1(trainingset2))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor',  'none', ...
    'LineStyle', '-', 'LineWidth', .5, ...
    'EdgeColor', PlottingParams.SeedColor,...
    'curvature', [.98 .1])

xlim([trainingNeurons{1}.candLat(1)-1 trainingNeurons{1}.candLat(end)+.5])
ylim([-7 9])
axis off;
set(gca, 'color', 'none')

%Rasters%
bottom = .1;
height = .45;
scale = .005;
spacing = .75/(2*nplots);

Red = trainingNeurons{1}.nIDs;
Green = trainingNeurons{2}.nIDs;
IsTrain1 = zeros(1,length(xdyn)); IsTrain1(Red) = 1;
IsTrain2 = zeros(1,length(xdyn)); IsTrain2(Green) = 1;

% keep track of which neurons participated in each syllable
FireDur1 = Latency{1}.FireDur(indkeep);
FireDur2= Latency{2}.FireDur(indkeep);

% classify neurons as specific or shared
Specific1 = FireDur1&~FireDur2;
Specific2 = FireDur2&~FireDur1;
Shared = (FireDur1&FireDur2);
indshared = find(Shared);

%collecting what I'll plot for the raster
sylIDtoplot = 3; % (don't choose a protosyllable
                 %  that's at the beginning of a bout)
k = length(union(trainingset1, trainingset2));
tindplot1 = trainingNeurons{1}.tind(sylIDtoplot) + ...
    trainingNeurons{1}.candLat-1; % time of example syl 1
tindplot2 = trainingNeurons{2}.tind(sylIDtoplot) + ...
    trainingNeurons{2}.candLat-1; % time of example syl 2
[~,indsort] = ...
    (sortrows(xdyn(:,[tindplot1 tindplot2]))); % sort by which fired first
tmp = xdyn(flipud(indsort), ...
    [tindplot1 tindplot2]); % pull out example from xdyn
IsTrain1 = IsTrain1(flipud(indsort));
IsTrain2 = IsTrain2(flipud(indsort));
Specific1 = Specific1(flipud(indsort));
Specific2 = Specific2(flipud(indsort));

indShared = (sum(tmp(:,1:length(tindplot1)),2)>0) & ...
    (sum(tmp(:,(length(tindplot1)+1):end),2)>0);
indBO = (sum(tmp(:,1:length(tindplot1)),2)>0) & ...
    (sum(tmp(:,(length(tindplot1)+1):end),2)==0);
rest = ~indShared;
tmp = tmp([find(indShared); find(rest)],:); % everything that will be 
                                            %  plotted in the rasters
IsTrain1 = IsTrain1([find(indShared); find(rest)]);
IsTrain2 = IsTrain2([find(indShared); find(rest)]);
Specific1 = Specific1([find(indShared); find(rest)]);
Specific2 = Specific2([find(indShared); find(rest)]);

% plot Bout Onset syllable
subplot('position', ...
    [ploti/nplots-2*spacing, bottom, length(tindplot1)*scale, height])
tmp1 = tmp(:,1:length(tindplot1)); % just bout onset syllable
PlottingParams.axesPosition = ...
    [ploti/nplots-2*spacing, bottom, length(tindplot1)*scale, height];

tOffset = trainingNeurons{1}.candLat(1)-1;
for j=1:size(tmp1,2) % for all the time steps
    Idx = find(tmp1(1:end-1,j)>0); % find the indices of active neurons
    if ~isempty(Idx)
        for k=1:length(Idx) % for all the active neurons
            Color = Specific1(Idx(k))*PlottingParams.Syl1Color + ...
                Specific2(Idx(k))*PlottingParams.Syl2Color;
            h = patch(10*([j-1,j,j,j-1]+tOffset),...
                [Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],...
                Color,'edgecolor','none');
        end
    end
end

hold on; box off
set(gca, 'fontsize', numFontSize)
set(gca, 'color', 'none', 'xtick', [0 100], ...
    'xticklabel', {'-100', '0'},'ydir', 'reverse', 'fontsize', numFontSize)
set(gca, 'ydir', 'reverse','tickdir','out','ticklength',[0.015 0.015],...
    'color', 'none', 'fontsize', numFontSize,'tickdir','out');

if PlottingParams.thisPanel==1
    ylabel('Neuron #','fontsize', labelFontSize)
    set(gca,'ytick',0:20:100,'fontsize', numFontSize)
else
    set(gca,'ytick',0:20:100,'yticklabel',{});
end

if PlottingParams.Hor
    if sum(indShared)>0 % if shared neurons
        plot([-10+trainingNeurons{1}.candLat(1)*10 ...
            trainingNeurons{1}.candLat(end)*10], ...
            (sum(indShared))*ones(1,2), ...
            'k', 'linewidth', PlottingParams.linewidth);
    end
    plot([-10+trainingNeurons{1}.candLat(1)*10 ...
        trainingNeurons{1}.candLat(end)*10], ...
        (sum(indBO)+sum(indShared))*ones(1,2), ...
        'k', 'linewidth', PlottingParams.linewidth);
end
if isfield(PlottingParams,'boutOnsetElement')
    patch([-10+trainingNeurons{1}.candLat(1)*10 -10 -10 ...
        -10+trainingNeurons{1}.candLat(1)*10],[-4 -4 -2 -2],Syl1Color);
    patch([-10+trainingNeurons{2}.candLat(1)*10 ...
        trainingNeurons{2}.candLat(end)*10 ...
        trainingNeurons{2}.candLat(end)*10 ...
        -10+trainingNeurons{2}.candLat(1)*10],[-4 -4 -2 -2],Syl2Color);
    if PlottingParams.thisPanel>1 
        text((-10+trainingNeurons{1}.candLat(1)*10-10)/2,-10,...
            '\epsilon','fontsize',7)
        text(((-10+trainingNeurons{2}.candLat(1)*10)+...
            (trainingNeurons{2}.candLat(end)*10))/2,-10,...
            '\alpha','fontsize',7)
    end
else
    patch([-10+trainingNeurons{1}.candLat(1)*10 ...
        trainingNeurons{1}.candLat(end)*10 ...
        trainingNeurons{1}.candLat(end)*10 ...
        -10+trainingNeurons{1}.candLat(1)*10],...
        [-4 -4 -2 -2],PlottingParams.Syl1BarColor);
    if PlottingParams.thisPanel>2
        text(((-10+trainingNeurons{1}.candLat(1)*10)+...
            (trainingNeurons{1}.candLat(end)*10))/2,-10,...
            '\beta','fontsize',7)
    end
end
ylim([-5 ntot])
xlim([-10+trainingNeurons{1}.candLat(1)*10 ...
    trainingNeurons{1}.candLat(end)*10+10])

% Plot protosyllable.
subplot('position', ...
    [ploti/nplots-spacing, bottom, length(tindplot2)*scale, height])
tmp1 = tmp(:,(length(tindplot1)+1):end); % just for protosyllable
PlottingParams.axesPosition = ...
    [ploti/nplots-spacing, bottom, length(tindplot2)*scale, height];

tOffset = trainingNeurons{2}.candLat(1)-1;
for j=1:size(tmp1,2) % for all the time steps
    Idx = find(tmp1(1:end-1,j)>0); % find the indices of active neurons
    if ~isempty(Idx)
        for k=1:length(Idx) % for all the active neurons
            Color = Specific1(Idx(k))*PlottingParams.Syl1Color + ...
                Specific2(Idx(k))*PlottingParams.Syl2Color;
            h = patch(10*([j-1,j,j,j-1]+tOffset),...
                [Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],...
                Color,'edgecolor','none');
        end
    end
end

hold on; box off
set(gca, 'fontsize', numFontSize)
set(gca, 'color', 'none', 'xtick', [0 100], ...
    'xticklabel', {'-100', '0'}, 'ydir', 'reverse', ...
    'fontsize', numFontSize)
set(gca, 'ydir', 'reverse','tickdir','out','ticklength',[0.015 0.015], ...
    'color', 'none', 'fontsize', numFontSize,'tickdir','out');
set(gca, 'ytick',0:20:100,'yticklabel',{})
if PlottingParams.Hor
    if sum(indShared)>0 % if shared neurons
        plot([-10+trainingNeurons{2}.candLat(1)*10 ...
            trainingNeurons{2}.candLat(end)*10], ...
            (sum(indShared))*ones(1,2), ...
            'k', 'linewidth', PlottingParams.linewidth);
    end
    plot([-10+trainingNeurons{2}.candLat(1)*10 ...
        trainingNeurons{2}.candLat(end)*10], ...
        (sum(indBO)+sum(indShared))*ones(1,2), ...
        'k', 'linewidth', PlottingParams.linewidth);
end

patch([-10+trainingNeurons{2}.candLat(1)*10 ...
    trainingNeurons{2}.candLat(end)*10 ...
    trainingNeurons{2}.candLat(end)*10 ...
    -10+trainingNeurons{2}.candLat(1)*10],...
    [-4 -4 -2 -2],PlottingParams.Syl2BarColor);
if isfield(PlottingParams,'boutOnsetElement')
    if PlottingParams.thisPanel>1
        text(((-10+trainingNeurons{2}.candLat(1)*10)+...
            (trainingNeurons{2}.candLat(end)*10))/2,-10,...
            '\alpha','fontsize',7)
    end
else
    if PlottingParams.thisPanel>2
        text(((-10+trainingNeurons{2}.candLat(1)*10)+...
            (trainingNeurons{2}.candLat(end)*10))/2,-10,...
            '\alpha','fontsize',7)
    end
end
ylim([-5 ntot])
xlim([-10+trainingNeurons{2}.candLat(1)*10 ...
    trainingNeurons{2}.candLat(end)*10+10])


end

%% Plotting functions: *plotSubsong*
function plotSubsong(w, xdyn, trainingNeurons, PlottingParams)
% Makes network diagram and raster plots, called by
% AlternatingDifferentiation
% w: weight matrix
% xdyn: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing 
%  neuron and time indices for each training neuron type
% PlottingParams: sets linewidth, etc.  

%Plotting parameters%
msize = PlottingParams.msize;
linewidth = PlottingParams.linewidth;
Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
ProtoSylColor = PlottingParams.ProtoSylColor;
numFontSize = PlottingParams.numFontSize;
labelFontSize = PlottingParams.labelFontSize;
nplots = PlottingParams.totalPanels; 
ploti = PlottingParams.thisPanel; 

%Network diagram%
subplot('position', [ploti/nplots-.7/nplots, .6, .7/nplots, .35])
cla; hold on

% calculate latency of each neuron
Latency = findLatency(xdyn, trainingNeurons);

% first keep track of latencies, and exclude neurons that don't fire at a
%  consistent latency
nsteps = size(xdyn,2); 
n = size(xdyn,1); 
ntot = n;
x = zeros(1,n);
y = zeros(1,n);
trainingset1 = trainingNeurons{1}.nIDs;
for ni = 1:n
    if length(intersect(trainingset1,ni))>0 % if it's a training neuron
        x(ni) = trainingNeurons{1}.candLat(1);
    else
        if Latency{1}.FireDur(ni) % if it participated in the syllable
            x(ni) = Latency{1}.mode(ni);
        else % if it didn't fire at consistent latency
            x(ni) = NaN;
        end
    end
end
indkeep = find(~isnan(x));
y = y(indkeep); 
w = w(indkeep,indkeep); 
xdyn = xdyn(indkeep,:); 
x = x(indkeep); 
ux = unique(x); 

% keep track of which neurons participated
FireDur1 = Latency{1}.FireDur(indkeep); 

% for each latency (x), spread along y
y1 = zeros(1,size(w,1)); 
for ui = 1:length(ux)
    indshared = (x==ux(ui))&FireDur1;
    [~,y1(indshared)] = sort(y1(indshared));
    tocentershared = 1+(numel(find(indshared))-1)/2;
    y1(indshared) = y1(indshared)-tocentershared;
end

% keep only feedforward part of weight matrix
wplot = w; 
n = size(wplot,1); 
for i = 1:n
    for j = 1:n
        ff = x(i)<x(j);
        longrange = abs(x(i)-x(j))>2; 
        if (~ff) | longrange
            wplot(j,i) = 0; 
        end
    end
end

% Color weights white to black between wplotmin and wplotmax
offset = .5;
wplot = wplot-PlottingParams.wplotmin+offset; 
wplot(wplot<0) = 0;
wplot = wplot/(PlottingParams.wplotmax-PlottingParams.wplotmin);
wplot(wplot<prctile(wplot(:), PlottingParams.wprctile)) = 0; 
wplotold = wplot; 
for i = 1:size(wplot,1)
    [~,ind] = sort(wplot(:,i), 'descend'); 
    indplot = zeros(size(wplot,1),1); 
    indplot(ind(1:min(PlottingParams.wperneuron,length(ind)))) = 1;   
    wplot(~indplot,i) = 0; 
end
for i = 1:size(wplot,1) 
    if sum(wplot(i,:)>0)<PlottingParams.wperneuronIn
        [~,ind] = sort(wplotold(i,:), 'descend'); 
        indplot = zeros(size(wplot,1),1); 
        wplot(i,ind(1:min(PlottingParams.wperneuron,length(ind)))) = ...
            wplotold(i,ind(1:min(PlottingParams.wperneuron,length(ind))));
    end
end

% jitter a little in x and y, so it doesn't look like a grid
%  but don't jitter seed neurons
jitter = .1; 
Seed0 = randn(1,500);
indJitter = setdiff(1:length(x), trainingset1); 
x(indJitter)= x(indJitter)+jitter*Seed0(1:length(x(indJitter)));
y1(indJitter) = y1(indJitter)+...
    jitter*Seed0((length(x(indJitter))+1):(2*length(x(indJitter))));


% plot w in order from weakest to strongest, so darker lines are on top
n = size(wplot,1); 
js = repmat((1:n)',1,n); 
is = repmat((1:n),n,1); 
isVec = is(:);
jsVec = js(:); 
wVec = wplot(:); 
[wSort,indSort] = sort(wVec, 'ascend'); 
for k = 1:length(wSort)
    i = isVec(indSort(k)); 
    j = jsVec(indSort(k)); 
    if wplot(j,i)>0
        ff = x(i)<=x(j);
        longrange = abs(x(i)-x(j))>2; 
        loopback = (round(x(i))==round(max(x)))&...
            (round(x(j))==round(min(x)));
        if ff & ~longrange%|loopback
            C = ones(1,3)-wplot(j,i)*ones(1,3);
            plot([x(i), x(j)], [y1(i),y1(j)], ...
                'color', C, 'linewidth', linewidth)
        end
    end
end

dotColor = zeros(length(x),3); 

for pli = 1:length(x)
    plot(x(pli),y1(pli), 'marker', '.', ...
        'color', dotColor(pli,:), 'markersize', msize)
end

% plot rectangle for seed neurons
rx = 1-.5;
ry = min(y1(trainingset1))-.5;
rw = 1;
rh = max(y1(trainingset1)) - min(y1(trainingset1))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor',  'none', ...
    'LineStyle', '-', 'LineWidth', .5, ...
    'EdgeColor',PlottingParams.SeedColor, ...
    'curvature', [.98 .1])

xlim([trainingNeurons{1}.candLat(1)-1 ...
    trainingNeurons{1}.candLat(end)/2+.5])
ylim([min(y1)-1 max(y1)+2])
axis off; 
set(gca, 'color', 'none')

%Rasters%
Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
ProtoSylColor = PlottingParams.ProtoSylColor;
numFontSize = PlottingParams.numFontSize;
labelFontSize = PlottingParams.labelFontSize;

bottom = .1; 
height = .45; 
scale = .005; 
spacing = .75/(2*nplots); 

%collecting what I'll plot for the raster
sylIDtoplot = 1; 
k = length(trainingset1); 
tindplot1 = trainingNeurons{1}.tind(sylIDtoplot) + ...
    trainingNeurons{1}.candLat-1; % time of example syl 1
[~,indsort] = (sortrows(xdyn(:,[tindplot1]))); % sort by which fired first
tmp = xdyn(flipud(indsort), [tindplot1]); % pull out example from xdyn

% plot raster 
subplot('position', ...
    [ploti/nplots-1.4*spacing, bottom, length(tindplot1)*scale, height])
tmp1 = tmp(:,1:length(tindplot1)); 
IsTrain = zeros(size(tmp,1)); IsTrain(trainingNeurons{1}.nIDs) = 1; 
tOffset = trainingNeurons{1}.candLat(1)-1; 
for j=1:size(tmp1,2) % for all the time steps
    Idx = find(tmp1(1:end-1,j)>0); % find the indices of active neurons    
    if ~isempty(Idx)
        for k=1:length(Idx) % for all the active neurons
            Color = IsTrain(Idx(k))*PlottingParams.SubsongSylColor;
            h = patch(10*([j-1,j,j,j-1]+tOffset),...
                [Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],...
                Color,'edgecolor','none');
        end  
    end
end

hold on; box off
set(gca, 'fontsize', numFontSize)
set(gca, 'color', 'none', 'xtick', 0:50:200,...
    'xticklabel',{'0','','100','','200'}, ...
    'ytick',0:20:100,'ydir', 'reverse', ...
    'tickdir','out','ticklength',[0.015 0.015],'fontsize', numFontSize)
ylabel('Neuron #','fontsize', labelFontSize)
ylim([-5 ntot])
xlim([0 trainingNeurons{1}.candLat(end)*10+10])
end