% Emily Mackevicius 12/10/2014, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

% plotting setup
clf;
clear all;

isEPS = 0; 

if isEPS 
    PlottingParams.msize = 8; % change to what is best for EPS figure
    PlottingParams.linewidth = .25; 
else
    PlottingParams.msize = 3;
    PlottingParams.linewidth = .25; 
end
PlottingParams.SeedColor = [1 .9 1];
PlottingParams.Syl2Color = [1 0 0]; 
PlottingParams.Syl1Color = [0 0 1]; 
PlottingParams.Syl2BarColor = [1 0 0]; 
PlottingParams.Syl1BarColor = [0 0 1];
PlottingParams.ProtoSylColor = [0 0 0]; 
PlottingParams.ProtoSylBarColor = [.5 .5 .5];
PlottingParams.numFontSize = 5; 
PlottingParams.labelFontSize = 8; 
PlottingParams.wplotmin = 0; 
PlottingParams.wplotmax = 2; % this should be wmaxSplit
PlottingParams.wprctile = 0; % plot all weights above this percentile. 
PlottingParams.totalPanels = 3; 
PlottingParams.thisPanel = 1; 
PlottingParams.sortby = 'weightMatrix'; 
PlottingParams.boutOnsetElement = 1; 

% Alternating seed neuron differentiation
figure(1); clf
set(gcf, 'color', ones(1,3));
if isEPS
    set(gcf, 'units','centimeters', 'position', [5 5 14 9])
end

seed = 3010; %978, 1009, 1012, 1021,1022, 1023
p.seed = seed; 
p.wmax = 1;             % single synapse hard bound
p.m = 5;                % desired number of synapses per neuron (wmax = Wmax/m)
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 500;         % time-steps to simulate -- each time-step is 1 burst duration.
p.pn = .01;             % probability of external stimulation of at least one neuron at any time
p.trainingInd = 1:10;   % index of training neurons
p.beta = .13;           % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .05;            % learning rate parameter
p.epsilon = .15;        % relative strength of heterosynaptic LTD
p.tau = 4;              % time constant of adaptation
p.gamma= .01;           % strength of recurrent inhibition
wmaxSplit = 2;          % single synapse hard bound to induce splitting (increased to encourage fewer stronger synapses)
gammaSplit =.05;        % increased strength of recurrent inhibition to induce splitting

Niter = [5    95   30   500]; % number of iterations for each plot (first 2 are protosyll, last 2 are splitting)
gammas = sigmf(1:Niter(end),[1/200 250])*gammaSplit; % gradually increase gamma to gammaSplit
p.gammas = gammas;
p.wmaxSplit = wmaxSplit; 
p.gammaSplit = gammaSplit; 
p.Niter = Niter; 

if ~isEPS
    folder = 'C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\SavedParams';
    timestamp = datestr(now, 'mmm-dd-yyyy-HH-MM-SS');
%     SavedHere = fullfile(folder, ['Params', timestamp])
%     save(SavedHere,'p');
end

PlotIters = 0; % set to 1, and increase Niter(3), if you want to plot each step as it goes

figure(1); clf

Wmax = p.wmax*p.m;
eta = p.eta; 

% random initial weights
rng(seed);
w0 = 2*rand(p.n)*Wmax/p.n; 
%
% training inputs
k = length(p.trainingInd);
trainint = p.trainint;
nsteps = p.nsteps;
n = p.n;
pn = p.pn;
% training inputs
CyclesPerBout = 5; 
bOnOffset = 6; 
HowClamped = 10; 
HowOn = 25; 
HowOnPsyl = 25; 
trainingNeurons{1}.nIDs = 1:k/2;
trainingNeurons{2}.nIDs = (k/2+1):k;
trainingNeurons{1}.candLat = (-bOnOffset+1):p.trainint;
trainingNeurons{2}.candLat =  1:p.trainint; 
trainingNeurons{1}.thres = 4;
trainingNeurons{2}.thres = 6;
Input = -HowClamped*ones(k, nsteps); % clamp training neurons
bOnOffsetVar = [1 randperm(20)];
indPsyl = [];
indBstart = [];
indOff = [];
prevPsylEnd = 1; 
for i = 1:(nsteps/CyclesPerBout/trainint)
    istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
    indPsyl = [indPsyl istart istart+trainint istart+2*trainint];
    indBstart = [indBstart istart-bOnOffset]; 
    indOff = [indOff prevPsylEnd:(istart-bOnOffset-1)];
    prevPsylEnd = istart+3*trainint;
end
indPsyl = indPsyl(indPsyl<=nsteps);
indBstart = indBstart(indBstart<=nsteps);
trainingNeurons{1}.tind = indBstart+bOnOffset;
trainingNeurons{2}.tind = setdiff(indPsyl, indBstart+bOnOffset); 
Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
Input(:,indOff) = -HowClamped; % clamp all neurons between bouts
trainingNeurons{1}.candLat = (-bOnOffset+1):trainint;
trainingNeurons{2}.candLat =  1:trainint; 
bdyn = double(rand(n,nsteps)>=(1-pn));
bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
bdyn(1:k,:) = Input; 
probeInput = bdyn; 
%%

w = w0; 
PlottingParams.thisPanel = 1;
niter = Niter(1);     % number of iterations to run
for j = 1:niter
    % Construct input
    Input = -HowClamped*ones(k, nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)];
    pSylVar = ceil(rand(1,20)*10);
    indPsyl = [];
    indBstart = [];
    indOff = [];
    prevPsylEnd = 1; 
    for i = 1:(nsteps/CyclesPerBout/trainint)
        istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
        indPsyl = [indPsyl [istart istart+trainint istart+2*trainint]];%+pSylVar(i)-bOnOffset];
        indBstart = [indBstart istart-bOnOffset]; 
        indOff = [indOff prevPsylEnd:(istart-bOnOffset-1)];
        prevPsylEnd = istart+3*trainint;
    end
    indPsyl = indPsyl(indPsyl<=nsteps);
    indBstart = indBstart(indBstart<=nsteps);
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    Input(:,indOff) = -HowClamped; % clamp all neurons between bouts
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    % One 'bout' of learning
    [w xdyn] = HVCIter(p);
end
%

p.eta = 0; p.input = probeInput; 
[w xdyn] = HVCIter(p); % probe run
p.eta = eta; 

PlottingParams.thisPanel = 1;
PlottingParams.Hor = 0; 
plotHVCnet_boutOnset(w, xdyn, trainingNeurons, PlottingParams)
PlottingParams.Hor = 1;

%% finish forming protosyllable
PlottingParams.thisPanel = 2;
niter = Niter(2);     % number of iterations to run
for j = 1:niter
    % Construct input
    Input = -HowClamped*ones(k, nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)];
    pSylVar = ceil(rand(1,20)*10);
    indPsyl = [];
    indBstart = [];
    indOff = [];
    prevPsylEnd = 1; 
    for i = 1:(nsteps/CyclesPerBout/trainint)
        istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
        indPsyl = [indPsyl [istart istart+trainint istart+2*trainint]];%+pSylVar(i)-bOnOffset];
        indBstart = [indBstart istart-bOnOffset]; 
        indOff = [indOff prevPsylEnd:(istart-bOnOffset-1)];
        prevPsylEnd = istart+3*trainint;
    end
    indPsyl = indPsyl(indPsyl<=nsteps);
    indBstart = indBstart(indBstart<=nsteps);
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    Input(:,indOff) = -HowClamped; % clamp all neurons between bouts
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    % One 'bout' of learning
    [w xdyn] = HVCIter(p);
end

p.eta = 0; p.input = probeInput; 
[w xdyn] = HVCIter(p); % probe run
p.eta = eta; 

PlottingParams.thisPanel = 2;
plotHVCnet_boutOnset(w, xdyn, trainingNeurons, PlottingParams)
wpsyl = w; 

%% splitting 
shg
w = wpsyl; 
p.wmax = wmaxSplit;  
p.m = Wmax/p.wmax;
PlottingParams.thisPanel = 3;
niter = Niter(3); 
for j = 1:niter
    % Construct input
    Input = -HowClamped*ones(k, nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)];
    pSylVar = ceil(rand(1,20)*10);
    indPsyl = [];
    indBstart = [];
    indOff = [];
    prevPsylEnd = 1; 
    for i = 1:(nsteps/CyclesPerBout/trainint)
        istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
        indPsyl = [indPsyl [istart istart+trainint istart+2*trainint]];%+pSylVar(i)-bOnOffset];
        indBstart = [indBstart istart-bOnOffset]; 
        indOff = [indOff prevPsylEnd:(istart-bOnOffset-1)];
        prevPsylEnd = istart+3*trainint;
    end
    indPsyl = indPsyl(indPsyl<=nsteps);
    indBstart = indBstart(indBstart<=nsteps);
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    Input(:,indOff) = -HowClamped; % clamp all neurons between bouts
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    p.gamma = gammas(j); 
    [w xdyn] = HVCIter(p);
    if  PlotIters & (mod(j,50)==0); % if you want to plot each step as it goes
        j
        subplot(1,4,3)
        plotHVCnet_boutOnset(w, xdyn, trainingNeurons, PlottingParams)
        pause(.5)
    end
end

p.eta = 0; p.input = probeInput; 
[w xdyn] = HVCIter(p); % probe run
p.eta = eta; 

% PlottingParams.thisPanel = 3;
% plotHVCnet_boutOnsetElement(w, xdyn, trainingNeurons, PlottingParams)
%%
PlottingParams.thisPanel = 4;
% Later splitting 
niter = Niter(4); 
for j = (Niter(3)+1):Niter(4)
    % Construct input
    Input = -HowClamped*ones(k, nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)];
    pSylVar = ceil(rand(1,20)*10);
    indPsyl = [];
    indBstart = [];
    indOff = [];
    prevPsylEnd = 1; 
    for i = 1:(nsteps/CyclesPerBout/trainint)
        istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
        indPsyl = [indPsyl [istart istart+trainint istart+2*trainint]];%+pSylVar(i)-bOnOffset];
        indBstart = [indBstart istart-bOnOffset]; 
        indOff = [indOff prevPsylEnd:(istart-bOnOffset-1)];
        prevPsylEnd = istart+3*trainint;
    end
    indPsyl = indPsyl(indPsyl<=nsteps);
    indBstart = indBstart(indBstart<=nsteps);
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    Input(:,indOff) = -HowClamped; % clamp all neurons between bouts
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(:,indOff) = -HowClamped; % clamp all neurons between bouts
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    p.gamma = gammas(j); 
    [w xdyn] = HVCIter(p);
end

p.eta = 0; p.input = probeInput; 
[w xdyn] = HVCIter(p); % probe run
p.eta = eta; 

PlottingParams.thisPanel = 3;
plotHVCnet_boutOnset(w, xdyn, trainingNeurons, PlottingParams)


%%
if isEPS
    cd('Z:\Fee_lab\Papers\HVC_differentiation\Figures\EPS_files');
    export_fig(1,'SuppFig10e.eps','-transparent','-eps','-painters');
else
    %figure parameters, exporting
    figw = 6*3/4;
    figh = 4; 
    set(gcf, 'color', [1 1 1],'papersize', [figw figh], 'paperposition', [0 0 figw*.9 figh])
    suptitle(['seed ', num2str(seed)])
    print -dmeta -r150
end