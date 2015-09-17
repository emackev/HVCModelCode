% Emily Mackevicius 12/10/2014, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

% plotting setup
clf;
clear;

isEPS = 0; 

if isEPS 
    PlottingParams.msize = 8; % change to what is best for EPS figure
    PlottingParams.linewidth = .25;
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultTextFontName', 'Arial')
    PlottingParams.labelFontSize = 7; 
else
    PlottingParams.msize = 5;
    PlottingParams.linewidth = 1;
    PlottingParams.labelFontSize = 6; 
end

Margin = 1/5; 
nplots = 4;
plotw = .23; 
netw = plotw-.01;
rasterw = plotw-Margin/2;
rasterh = 3/4; 
netoffset = Margin/3;
neth = 1/4-Margin/4-.01; 
PlottingParams.Syl1Color = [1 0 0]; 
PlottingParams.Syl2Color = [0 1 0]; % please choose orthogonal colors.. if you don't I'll try and normalize colors and it'll look muddy
PlottingParams.ProtoSylColor = [1 0 1]; 
PlottingParams.Syl1Color = PlottingParams.Syl1Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
PlottingParams.Syl2Color = PlottingParams.Syl2Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
PlottingParams.numFontSize = 5; 
PlottingParams.wplotmin = 0; 
PlottingParams.wplotmax = 2; % this should be wmaxSplit
PlottingParams.totalPanels = 4; 
PlottingParams.sortby = 'activity'; 

PlotIters = 0; 
% Alternating seed neuron differentiation

seed = 210
p.seed = seed;          % seed random number generator
p.wmax = 1;             % single synapse hard bound
p.m = 9;               % desired number of synapses per neuron (wmax = Wmax/m)
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 100;         % time-steps to simulate -- each time-step is 1 burst duration.
p.pn = .01;             % probability of external stimulation of at least one neuron at any time
p.trainingInd = 1:9;   % index of training neurons
p.beta = .115;          % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .025;           % learning rate parameter
p.epsilon = .2;         % relative strength of heterosynaptic LTD
p.tau = 4;              % time constant of adaptation
p.gamma= .01;           % strength of recurrent inhibition

wmaxSplit = 3;          % single synapse hard bound to induce splitting (increased to encourage fewer stronger synapses)
gammaSplit =.18;        % increased strength of recurrent inhibition to induce splitting


Niter = [1 499 690 2000]; % number of iterations for each plot (first 2 are protosyll, last 2 are splitting)
gammas = sigmf(1:Niter(end),[1/200 500])*gammaSplit; % gradually increase gamma to gammaSplit
Wmax = p.wmax*p.m;

% saving params for later.
p.gammas = gammas;
p.wmaxSplit = wmaxSplit; 
p.gammaSplit = gammaSplit; 
p.Niter = Niter; 

if ~isEPS
    folder = 'C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\SavedParams';
    timestamp = datestr(now, 'mmm-dd-yyyy-HH-MM-SS');
    SavedHere = fullfile(folder, ['Params', timestamp])
    save(SavedHere,'p');
end
%%
% random initial weights
rng(seed);
w0 = 2*rand(p.n)*Wmax/p.n; 

%Psyl inputs
% training inputs
k = length(p.trainingInd);
trainint = p.trainint;
nsteps = p.nsteps;
n = p.n;
pn = p.pn;
HowClamped = 10; 
HowOn = 10; 
Input = -HowClamped*ones(k, nsteps); %clamp training neurons (effectively giving them higher threshold)
Input(:,mod(1:nsteps,trainint)==1) = HowOn; % rhythmic activation of training neurons
PsylInput = Input; 

%Alternating Inputs
% training inputs
Input =-HowClamped*ones(k, nsteps); % clamp training neurons (effectively giving them higher threshold)
Input(1:3,mod(1:nsteps,3*trainint)==1) = HowOn; % alternating rhythmic activation of training neurons
Input(4:6,mod(1:nsteps,3*trainint)==trainint+1) = HowOn; % alternating rhythmic activation of training neurons
Input(7:9,mod(1:nsteps,3*trainint)==2*trainint+1) = HowOn; % alternating rhythmic activation of training neurons
AltInput = Input;

figure(1); clf
set(gcf, 'color', ones(1,3));
if isEPS
    set(gcf, 'units','centimeters', 'position', [5 5 18 9])
end

PlottingParams.thisPanel = 1; 
w = w0;
niter = Niter(1);     % number of iterations to run
for i = 1:niter
    % Construct input
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = PsylInput; 
    p.w = w; 
    p.input = bdyn;
    % One 'bout' of learning
    %tmp = p; tmp.eta = 0; 
    [w xdyn] = HVCIter(p);
end
HVCtestRaster_intoThree(xdyn,PsylInput,w,PlottingParams);
wpsyl = w; 

PlottingParams.thisPanel = 2;

niter = Niter(2);     % number of iterations to run
for i = 1:niter
    % Construct input
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = PsylInput; 
    p.w = w; 
    p.input = bdyn;
    % One 'bout' of learning
    %tmp = p; tmp.eta = 0; 
    [w xdyn] = HVCIter(p);
end
HVCtestRaster_intoThree(xdyn,PsylInput,w,PlottingParams);
wpsyl = w; 

%%
PlottingParams.thisPanel = 3;

%  splitting 
w = wpsyl;
p.wmax = wmaxSplit;  
p.m = Wmax/p.wmax;

niter = Niter(3); 
for i = 1:niter
    % Construct input
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = AltInput; 
    p.w = w; 
    p.input = bdyn;
    p.gamma = gammas(i); 
    [w xdyn] = HVCIter(p);
    if  PlotIters &   mod(i,50)==0 ; % if you want to plot each step as it goes
        i
        HVCtestRaster(xdyn,AltInput,w,PlottingParams);
        %subplot('position', [netoffset+2*plotw Margin/4+rasterh netw neth]); plotHVCnet(w,xdyn,trainint,trainingNeurons,PlottingParams)
        %plotHVCnet(w, xdyn, trainint, trainingNeurons, PlottingParams);
        pause(.5)
    end
end
HVCtestRaster_intoThree(xdyn,AltInput,w,PlottingParams);

%%
PlottingParams.thisPanel = 4;
% Later splitting 
niter = Niter(4); 
for i = (Niter(3)+1):Niter(4)
    % Construct input
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = AltInput; 
    % One 'bout' of learning
    p.w = w; 
    p.input = bdyn;
    p.gamma = gammas(i); 
    [w xdyn] = HVCIter(p);
end
HVCtestRaster_intoThree(xdyn,AltInput,w,PlottingParams);


%%
if isEPS
    cd('Z:\Fee_lab\Papers\HVC_differentiation\Figures\EPS_files');
    export_fig(1,'Fig7a.eps','-transparent','-eps','-painters');
else
    figure parameters, exporting
    figw = 6;
    figh = 2;
    set(gcf, 'color', [1 1 1],'papersize', [figw figh], 'paperposition', [0 0 figw figh])
    suptitle(['seed ', num2str(seed)])
    print -dmeta -r150
end