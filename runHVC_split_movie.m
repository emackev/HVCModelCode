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

highQual = 0; 

% PlottingParams.msize = 25;
% PlottingParams.linewidth = 1; 
% PlottingParams.Syl1Color = [1 0 0]; 
% PlottingParams.Syl2Color = [0 1 0]; % please choose orthogonal colors.. if you don't I'll try and normalize colors and it'll look muddy
% PlottingParams.ProtoSylColor = [1 0 1]; 
% PlottingParams.Syl1Color = PlottingParams.Syl1Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
% PlottingParams.Syl2Color = PlottingParams.Syl2Color/max(PlottingParams.Syl1Color+PlottingParams.Syl2Color);
% PlottingParams.numFontSize = 5; 
% PlottingParams.labelFontSize = 14; 
% PlottingParams.wplotmin = 0; 
% PlottingParams.wplotmax = 2; % this should be wmaxSplit
% PlottingParams.wprctile = 0; % plot all weights above this percentile.  If nonzero, ignores wplotmin, wplotmax
% PlottingParams.wperneuron = 6; % max outgoing weights plotted
% PlottingParams.wperneuronIn = 9; % min incoming weights plotted
rng('default')

% Alternating seed neuron differentiation

seed = 9038;
p.seed = seed;          % seed random number generator
p.n = 100;              % n neurons
p.trainint = 10;        % Time interval between inputs
p.nsteps = 100;         % time-steps to simulate -- each time-step is 1 burst duration.
p.pn = .01;             % probability of external stimulation of at least one neuron at any time
p.trainingInd = 1:10;   % index of training neurons
p.beta = .115;          % strength of feedforward inhibition
p.alpha = 30;           % strength of neural adaptation
p.eta = .025;           % learning rate parameter
p.epsilon = .2;         % relative strength of heterosynaptic LTD
p.tau = 4;              % time constant of adaptation
gammaStart= .01;        % strength of recurrent inhibition
gammaSplit =.18;        % increased strength of recurrent inhibition to induce splitting
wmaxStart = 1;          % single synapse hard bound
wmaxSplit = 2;          % single synapse hard bound to induce splitting (increased to encourage fewer stronger synapses)
mStart = 10;            % desired number of synapses per neuron (wmax = Wmax/m)
Wmax = mStart*wmaxStart;% soft bound for weights of each neuron
mSplit = Wmax/wmaxSplit;% keep Wmax constant, change m & wmax to induce fewer stronger synapses
HowClamped = 10;        % give training neurons higher threshold
HowOn = 10;             % higher inputs to training neurons

nIterProto = 500;       % end of protosyllable stage
nIterPlotSplit1 = 492;  % number of splitting iterations before plotting intermediate splitting phase
nIterPlotSplit2 = 2000; % total number of splitting iterations

% parameters that change over development
protosyllableStage = [ones(1,nIterProto) zeros(1,nIterPlotSplit2)]; 
splittingStage = [zeros(1,nIterProto) ones(1,nIterPlotSplit2)];
gammas = gammaStart*ones(1,nIterProto);     % keep gamma at gammaStart during protosyllable stage
gammas = [gammas sigmf(1:nIterPlotSplit2,[1/200 500])*gammaSplit];  % gradually increase gamma to gammaSplit during splitting
wmaxs = protosyllableStage*wmaxStart + splittingStage*wmaxSplit;  
ms = protosyllableStage*mStart + splittingStage*mSplit;  

% plotting params and inputs 
nplots = 4;

PlottingParams.msize = 50; %20
PlottingParams.linewidth = 1; %.5;
PlottingParams.labelFontSize = 20; %12; 
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
PlottingParams.wprctile = 0; % plot all weights above this percentile.  If nonzero, ignores wplotmin, wplotmax
PlottingParams.wperneuron = 6; % max outgoing weights plotted
PlottingParams.wperneuronIn = 9; % min incoming weights plotted
bottom = .1; 
height = .55; 
scale = .005; 
spacing = .75/(2*nplots); 

%Psyl inputs
k = length(p.trainingInd);
Input = -HowClamped*ones(k, p.nsteps); %clamp training neurons (effectively giving them higher threshold)
Input(:,mod(1:p.nsteps,p.trainint)==1) = HowOn; % rhythmic activation of training neurons
PsylInput = Input; 
%Alternating Inputs
Input =-HowClamped*ones(k, p.nsteps); % clamp training neurons (effectively giving them higher threshold)
Input(1:k/2,mod(1:p.nsteps,2*p.trainint)==1) = HowOn; % alternating rhythmic activation of training neurons
Input((k/2+1):k,mod(1:p.nsteps,2*p.trainint)==p.trainint+1) = HowOn; % alternating rhythmic activation of training neurons
AltInput = Input;
%Subsong Inputs
nstepsSubsong = 1000; 
rng(seed)
isOnset = rand(1,nstepsSubsong)>.9; 
Input =-HowClamped*ones(k, nstepsSubsong); % clamp training neurons (effectively giving them higher threshold)
Input(1:k,isOnset) = HowOn; 
bdyn = double(rand(p.n,nstepsSubsong)>=(1-p.pn)); % Random activation
bdyn(1:k,:) = Input; 
subsongInput = bdyn; 

% plotting subsong 
for i = 1:2
    trainingNeurons{i}.nIDs = 1:k;
    trainingNeurons{i}.tind = find(isOnset);
    trainingNeurons{i}.candLat = 1:2*p.trainint; 
    trainingNeurons{i}.thres = 12; % criteria for participation during subsong (thres from testLatSig -- must fire at consistent latency more than 12 times in the bout of ~100 syllables to count as participating)
end
trainingNeuronsSubsong = trainingNeurons; clear trainingNeurons;

% plotting protosyl
trainingNeurons{1}.nIDs = 1:k;
trainingNeurons{2}.nIDs = 1:k;
trainingNeurons{1}.tind = find(mod(1:p.nsteps, p.trainint)==1);% repmat([true(1,p.trainint) false(1,p.trainint)],1,p.nsteps/p.trainint/2);
trainingNeurons{2}.tind = find(mod(1:p.nsteps, p.trainint)==1); %repmat([true(1,p.trainint) false(1,p.trainint)],1,p.nsteps/p.trainint/2);
trainingNeurons{1}.candLat = 1:p.trainint; 
trainingNeurons{2}.candLat = 1:p.trainint; 
trainingNeurons{1}.thres = 4; 
trainingNeurons{2}.thres = 4; 
trainingNeuronsPsyl = trainingNeurons; clear trainingNeurons;

% plotting splitting stages
trainingNeurons{1}.nIDs = 1:k/2;
trainingNeurons{2}.nIDs = (k/2+1):k;
trainingNeurons{1}.tind = find(mod(1:p.nsteps, 2*p.trainint)==1);% repmat([true(1,p.trainint) false(1,p.trainint)],1,p.nsteps/p.trainint/2);
trainingNeurons{2}.tind = find(mod(1:p.nsteps, 2*p.trainint)==p.trainint+1); %repmat([true(1,p.trainint) false(1,p.trainint)],1,p.nsteps/p.trainint/2);
trainingNeurons{1}.candLat = 1:p.trainint; 
trainingNeurons{2}.candLat = 1:p.trainint; 
trainingNeurons{1}.thres = 2; 
trainingNeurons{2}.thres = 2; 
trainingNeuronsAlt = trainingNeurons; clear trainingNeurons;


% set up to record movie
folder = fileparts(mfilename('fullpath')); %'C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\NetworkMovies';
timestamp = datestr(now, 'mmm-dd-yyyy-HH-MM-SS');
filename = ['NetLearnsSeed' num2str(seed) timestamp];
%aviobj = avifile(fullfile(folder, filename), 'compression', 'none', 'fps',
%20); -- old matlab version
obj = vision.VideoFileWriter(fullfile(folder, [filename, '.mp4']));
obj.FrameRate = 20; 
obj.FileFormat='MPEG4';
%obj.VideoCompressor='None (uncompressed)';
%% run simulation

% initialize weight matrix
rng(seed);
w0 = 2*rand(p.n)*Wmax/p.n; 

% subsong stage
pSubsong = p; 
pSubsong.gamma = gammas(1); 
pSubsong.wmax = wmaxs(1); 
pSubsong.m = ms(1); 
pSubsong.eta = 0; 
pSubsong.nsteps = nstepsSubsong; 
pSubsong.w = w0; 
pSubsong.input = subsongInput;
% Run subsong network
[wSubsong xdynSubsong] = HVCIter(pSubsong);
w = wSubsong;

figure
set(gcf, 'color', ones(1,3), 'units', 'inches', 'position', [.1 1 16 9]);
plotHVCnet(wSubsong, xdynSubsong, p.trainint, trainingNeuronsSubsong, PlottingParams)
text(5,10,'Subsong', 'fontsize', PlottingParams.labelFontSize, 'color', [0 0 0], 'horizontalalignment', 'center', 'verticalalignment', 'bottom')%, 'fontweight', 'bold')

set(gca, 'color', 'none');
if highQual
    s = rng; 
    %myaa % anti-aliasing
    rng(s);
end
F = getframe(gcf)
slowrate =20;
for l = 1:slowrate
    step(obj, F.cdata);
    %aviobj = addframe(aviobj,F); -- old matlab version
end
close all
%%
% learning stages
wOld = w; 
for i = 1:(nIterProto+nIterPlotSplit2)
    p.w = w;
    % set parameters that change over development
    p.gamma = gammas(i); 
    p.wmax = wmaxs(i); 
    p.m = ms(i); 
    % Construct input
    bdyn = double(rand(p.n,p.nsteps)>=(1-p.pn)); % Random activation
    bdyn(1:k,:) = protosyllableStage(i)*PsylInput+splittingStage(i)*AltInput; % drive seed neurons
    p.input = bdyn;
    % run one iteration
    [w xdyn] = HVCIter(p);
    dw(i) = norm(w(:) - wOld(:));
    wOld = w; 
    % add frames to movie
    if dw(i)>.1 | i == (nIterProto+nIterPlotSplit2)
        if i < 5
            slowrate = 20;
        else
            slowrate = 1;
        end
        figure
        set(gcf, 'color', ones(1,3), 'units', 'inches', 'position', [.1 1 16 9]);
        if protosyllableStage(i)
            plotHVCnet(w,xdyn,p.trainint,trainingNeuronsPsyl,PlottingParams)
            text(5,10,['Protosyllable stage: iteration ', num2str(i)], 'fontsize', PlottingParams.labelFontSize, 'color', [0 0 0], 'horizontalalignment', 'center', 'verticalalignment', 'bottom')%, 'fontweight', 'bold')
%             title(['Protosyllable stage: iteration ', num2str(i)], ...
%                 'fontsize', PlottingParams.labelFontSize, 'color', [1 0 1])%, 'fontweight', 'bold')
        else
            plotHVCnet(w,xdyn,p.trainint,trainingNeuronsAlt,PlottingParams)
            text(5,10,['Splitting stage: iteration ', num2str(i - nIterProto)], 'fontsize', PlottingParams.labelFontSize, 'color', [1 0 0], 'horizontalalignment', 'center', 'verticalalignment', 'bottom')%, 'fontweight', 'bold')
        end
        set(gca, 'color', 'none');
        if highQual
            s = rng; 
            myaa % anti-aliasing
            rng(s);
        end
        F = getframe(gcf) 
        for l = 1:slowrate
            step(obj, F.cdata);
            %aviobj = addframe(aviobj,F); -- old matlab version
        end
        close all
    end

%     
%     if i<5
%         slowrate = 20;
%         figure
%         set(gcf, 'color', ones(1,3));
%         plotHVCnet(w,xdyn,p.trainint,trainingNeuronsPsyl,PlottingParams)
%         title(['Protosyllable stage: iteration ', num2str(i)], 'fontsize', PlottingParams.labelFontSize)
%         set(gca, 'color', 'none');
%         if highQual
%             s = rng; 
%             myaa % anti-aliasing
%             rng(s);
%         end
%         F = getframe(gcf) 
%         for l = 1:slowrate
%             aviobj = addframe(aviobj,F);
%         end
%         close all
%     elseif i < nIterProto
%         speedrate = 10; 
%         if mod(i,speedrate)==0
%             figure
%             set(gcf, 'color', ones(1,3));
%             plotHVCnet(w,xdyn,p.trainint,trainingNeuronsPsyl,PlottingParams)
%             title(['Protosyllable stage: iteration ', num2str(i)], 'fontsize',PlottingParams.labelFontSize)
%             set(gca, 'color', 'none');
%             if highQual
%                 s = rng; 
%                 myaa % anti-aliasing
%                 rng(s);
%             end
%             F = getframe(gcf)
%             aviobj = addframe(aviobj,F);
%             close all
%         end
%     elseif i < 990
%         speedrate = 10; 
%         if mod(i,speedrate)==0
%             figure
%             set(gcf, 'color', ones(1,3));
%             plotHVCnet(w,xdyn,p.trainint,trainingNeuronsAlt,PlottingParams)
%             title(['Splitting stage: iteration ', num2str(i - nIterProto)], 'fontsize',PlottingParams.labelFontSize)
%             set(gca, 'color', 'none');
%             if highQual
%                 s = rng; 
%                 myaa % anti-aliasing
%                 rng(s);
%             end
%             F = getframe(gcf)
%             aviobj = addframe(aviobj,F);
%             close all
%         end
%     elseif i < 998 
%         slowrate = 20;
%         figure
%         set(gcf, 'color', ones(1,3));
%         plotHVCnet(w,xdyn,p.trainint,trainingNeuronsAlt,PlottingParams)
%         title(['Splitting stage: iteration ', num2str(i - nIterProto)], 'fontsize',PlottingParams.labelFontSize)
%         set(gca, 'color', 'none');
%         if highQual
%             s = rng; 
%             myaa % anti-aliasing
%             rng(s);
%         end
%         F = getframe(gcf) 
%         for l = 1:slowrate
%             aviobj = addframe(aviobj,F);
%         end
%         close all
%     elseif i < 1200
%         speedrate = 2; 
%         if mod(i,speedrate)==0
%             figure
%             set(gcf, 'color', ones(1,3));
%             plotHVCnet(w,xdyn,p.trainint,trainingNeuronsAlt,PlottingParams)
%             title(['Splitting stage: iteration ', num2str(i - nIterProto)], 'fontsize',PlottingParams.labelFontSize)
%             set(gca, 'color', 'none');
%             if highQual
%                 s = rng; 
%                 myaa % anti-aliasing
%                 rng(s);
%             end
%             F = getframe(gcf)
%             aviobj = addframe(aviobj,F);
%             close all
%         end
%     elseif i <= nIterProto + nIterPlotSplit2
%         speedrate = 20; 
%         if mod(i,speedrate)==0
%             figure
%             set(gcf, 'color', ones(1,3));
%             plotHVCnet(w,xdyn,p.trainint,trainingNeuronsAlt,PlottingParams)
%             title(['Splitting stage: iteration ', num2str(i - nIterProto)], 'fontsize',PlottingParams.labelFontSize)
%             set(gca, 'color', 'none');
%             if highQual
%                 s = rng; 
%                 myaa % anti-aliasing
%                 rng(s);
%             end
%             F = getframe(gcf)
%             aviobj = addframe(aviobj,F);
%             close all
%         end
%     end
end
slowrate = 20; % leave last network state on screen longer
for l = 1:slowrate
    step(obj, F.cdata);
    %aviobj = addframe(aviobj,F); -- old matlab version
end
release(obj);
