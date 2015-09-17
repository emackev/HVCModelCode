% Emily Mackevicius 12/10/2014, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

% plotting setup
clf;
clear all;
Margin = 1/5; 
nplots = 4;
plotw = .23; 
netw = plotw-.01;
rasterw = plotw-Margin/2;
rasterh = 3/4; 
netoffset = Margin/3;
neth = 1/4-Margin/4-.01; 
PlottingParams.msize = 5;
PlottingParams.linewidth = .01; 
PlottingParams.Syl1Color = [1 0 0]; 
PlottingParams.Syl2Color = [0 1 0]; % please choose orthogonal colors.. if you don't I'll try and normalize colors and it'll look muddy
PlottingParams.ProtoSylColor = [1 0 1]; 
PlottingParams.Syl1Color = PlottingParams.Syl1Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
PlottingParams.Syl2Color = PlottingParams.Syl2Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
PlottingParams.numFontSize = 5; 
PlottingParams.labelFontSize = 8; 
PlottingParams.wplotmin = 0; 
PlottingParams.wplotmax = 2; % this should be wmaxSplit
PlottingParams.totalPanels = 1; 
PlottingParams.thisPanel = 1; 

% Alternating seed neuron differentiation
figure(1); clf
set(gcf, 'color', ones(1,3));

seed = 8062
p.seed = seed;          % seed random number generator
p.wmax = 1;             % single synapse hard bound
p.m = 5;               % desired number of synapses per neuron (wmax = Wmax/m)
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 100;         % time-steps to simulate -- each time-step is 1 burst duration.
p.pn = .01;             % probability of external stimulation of at least one neuron at any time
p.trainingInd = 1:10;   % index of training neurons
p.beta = .1;           % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .025;           % learning rate parameter
p.epsilon = .1;         % relative strength of heterosynaptic LTD
p.tau = 3;              % time constant of adaptation
p.gamma= .01;           % strength of recurrent inhibition
wmaxSplit = 2;          % single synapse hard bound to induce splitting (increased to encourage fewer stronger synapses)
gammaSplit =.15;        % increased strength of recurrent inhibition to induce splitting

Niter = [1 1500 1500 1500]; % number of iterations for each plot (first 2 are protosyll, last 2 are splitting)
gammas = sigmf(1:Niter(end),[1/200 500])*gammaSplit; % gradually increase gamma to gammaSplit
p.gammas = gammas;
p.wmaxSplit = wmaxSplit; 
p.gammaSplit = gammaSplit; 
p.Niter = Niter; 



PlotIters = 1; % set to 1, and increase Niter(3), if you want to plot each step as it goes

figure(1); cla
% set up to record movie
folder = fileparts(mfilename('fullpath'));% 'C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\NetworkMovies';
timestamp = datestr(now, 'mmm-dd-yyyy-HH-MM-SS');
filename = ['NetLearnsSeed' num2str(seed) timestamp];
writerobj = VideoWriter(fullfile(folder, filename));
writerobj.FrameRate = 5; 
open(writerobj);


Wmax = p.wmax*p.m;

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
bOnOffset = 5; 
CyclesPerBout = 4; 
trainingNeurons{1}.nIDs = 1:k/2;
trainingNeurons{2}.nIDs = (k/2+1):k;
HowClamped = 10; 
HowOn = 10; 
Input = -HowClamped*ones(k, nsteps); %clamp training neurons (effectively giving them higher threshold)
bOnOffsetVar = randperm(20);  % decouple bout onset and protosyllables
indPsyl = [];
indBstart = []; 
Input(trainingNeurons{2}.nIDs,:) = -HowClamped;
for i = 1:(nsteps/CyclesPerBout/trainint)
    istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i); 
    indPsyl = [indPsyl istart istart+trainint istart+2*trainint];
    indBstart = [indBstart (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(end-i)]; 
end
Input(trainingNeurons{2}.nIDs,indPsyl) = HowOn; % alternating rhythmic activation of training neurons
Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons

imagesc(Input)

w = w0; 
%
% finish forming protosyllable
niter = Niter(2);     % number of iterations to run
for i = 1:niter
    % Construct input
    Input = -HowClamped*ones(k, nsteps); %clamp training neurons (effectively giving them higher threshold)
    bOnOffsetVar = randperm(20);  % decouple bout onset and protosyllables
    indPsyl = [];
    indBstart = []; 
    Input(trainingNeurons{2}.nIDs,:) = -HowClamped;
    for j = 1:(nsteps/CyclesPerBout/trainint)
        istart = (j-1)*CyclesPerBout*trainint+1+bOnOffsetVar(j); 
        indPsyl = [indPsyl istart istart+trainint istart+2*trainint];
        indBstart = [indBstart (j-1)*CyclesPerBout*trainint+1+bOnOffsetVar(end-j)]; 
    end
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOn; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    % One 'bout' of learning
    [w xdyn] = HVCIter(p);
    if (mod(i,50) == 1)|i<20
        HVCtestRaster_forMovies(xdyn,Input,w,PlottingParams)
        pause(.5)
        i
        frame = getframe(gcf);
        slowrate = 1; 
        for l = 1:slowrate
            writeVideo(writerobj,frame);
        end
        pause(.2)
    end
end
HVCtestRaster_forMovies(xdyn,Input,w,PlottingParams)
wpsyl = w; 


%%
% splitting 

shg
w = wpsyl; 
p.wmax = wmaxSplit;  
p.m = Wmax/p.wmax;

% training inputs

bOnOffset = 3; 
HowClamped = 10; 
HowOn = 10; 
HowOnPsyl = 1; 
trainingNeurons{1}.nIDs = 1:k/2;
trainingNeurons{2}.nIDs = (k/2+1):k;
Input = -HowClamped*ones(k, nsteps); % clamp training neurons
bOnOffsetVar = [1 randperm(20)];
indPsyl = [];
indBstart = [];
for i = 1:(nsteps/CyclesPerBout/trainint)
    istart = (i-1)*CyclesPerBout*trainint+1+bOnOffsetVar(i)+bOnOffset; 
    indPsyl = [indPsyl istart istart+trainint istart+2*trainint];
    indBstart = [indBstart istart-bOnOffset]; 
end
Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
imagesc(Input)

%

niter = Niter(3); 
for i = 1:niter
    % Construct input
    Input = -HowClamped*ones(k, nsteps); % clamp training neurons
    bOnOffsetVar = [1 randperm(20)];
    indPsyl = [];
    indBstart = [];
    for j = 1:(nsteps/CyclesPerBout/trainint)
        istart = (j-1)*CyclesPerBout*trainint+1+bOnOffsetVar(j)+bOnOffset; 
        indPsyl = [indPsyl istart istart+trainint istart+2*trainint];
        indBstart = [indBstart istart-bOnOffset]; 
    end
    Input(trainingNeurons{2}.nIDs,indPsyl) = HowOnPsyl; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{1}.nIDs,indBstart) = HowOn; % alternating rhythmic activation of training neurons
    
    bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
    bdyn(1:k,:) = Input; 
    p.w = w; 
    p.input = bdyn;
    p.gamma = gammas(i); 
    [w xdyn] = HVCIter(p);
%     Latency = findHVClatency(xdyn,trainint,trainingNeurons); 
%     Nsplit = sum(xor(Latency{1}.FireDur,Latency{2}.FireDur));
    if  PlotIters & (mod(i,10)==0); % if you want to plot each step as it goes
        i
        HVCtestRaster_forMovies(xdyn,Input,w,PlottingParams)
        frame = getframe(gcf);
        for l = 1:slowrate
            writeVideo(writerobj,frame);
        end
        pause(.2)
    end
end

cla;
HVCtestRaster_forMovies(xdyn,Input,w,PlottingParams)

% % Later splitting 
% niter = Niter(4); 
% for i = (Niter(3)+1):Niter(4)
%     % Construct input
%     bdyn = (rand(n,nsteps)>=(1-pn)); % Random activation
%     bdyn(1:k,:) = Input; 
%     % One 'bout' of learning
%     p.w = w; 
%     p.input = bdyn;
%     p.gamma = gammas(i); 
%     [w xdyn] = HVCIter(p);
% end
% 
% 
% HVCtestRaster(xdyn,Input)
% 
% 
% % figure parameters
% figw = 6;
% figh = 3; 
% set(gcf, 'color', [1 1 1],'papersize', [figw figh], 'paperposition', [0 0 figw*.9 figh])
% suptitle(['seed ', num2str(seed)])
% print -dmeta -r150

close(writerobj);
