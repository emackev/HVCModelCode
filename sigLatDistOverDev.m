clear all; close all
figure(1); clf
set(gcf, 'color', ones(1,3));


p.wmax = 1;             % single synapse hard bound
p.m = 10;               % desired number of synapses per neuron (wmax = Wmax/m)
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
p.gamma= .01;           % strength of recurrent inhibition
wmaxSplit = 2;          % single synapse hard bound to induce splitting (increased to encourage fewer stronger synapses)
gammaSplit =.18;        % increased strength of recurrent inhibition to induce splitting
Niter = [1 500 500]; 
gammas = sigmf(1:Niter(end),[1/200 500])*gammaSplit; % gradually increase gamma to gammaSplit

Wmax = p.wmax*p.m;
wmax = p.wmax; 
gamma = p.gamma; 
m = p.m; 

k = length(p.trainingInd);
trainint = p.trainint;
nsteps = p.nsteps;
n = p.n;
pn = p.pn;
HowClamped = 10; 
HowOn = 10; 

%Psyl inputs
% training inputs
trainingNeurons{1}.nIDs = 1:k;
trainingNeurons{2}.nIDs = 1:k;
trainingNeurons{1}.tind = find(mod(1:nsteps,trainint)==1); %
trainingNeurons{2}.tind = find(mod(1:nsteps,trainint)==1); %
trainingNeurons{1}.candLat = 1:trainint; 
trainingNeurons{2}.candLat = 1:trainint; 
trainingNeurons{1}.thres = 4; 
trainingNeurons{2}.thres = 4;
Input = -HowClamped*ones(k, nsteps); %clamp training neurons (effectively giving them higher threshold)
Input(:,mod(1:nsteps,trainint)==1) = HowOn; % rhythmic activation of training neurons
PsylInput = Input; 
trainingNeuronsPsyl = trainingNeurons; clear trainingNeurons;

%Alternating Inputs
% training inputs
trainingNeurons{1}.nIDs = 1:k/2;
trainingNeurons{2}.nIDs = (k/2+1):k;
trainingNeurons{1}.tind = find(mod(1:nsteps,2*trainint)==1); %repmat([true(1,trainint) false(1,trainint)],1,nsteps/trainint/2);
trainingNeurons{2}.tind = find(mod(1:nsteps,2*trainint)==trainint+1); %repmat([false(1,trainint) true(1,trainint)],1,nsteps/trainint/2);
trainingNeurons{1}.candLat = 1:trainint; 
trainingNeurons{2}.candLat = 1:trainint; 
trainingNeurons{1}.thres = 2; 
trainingNeurons{2}.thres = 2;
Input =-HowClamped*ones(k, nsteps); % clamp training neurons (effectively giving them higher threshold)
Input(trainingNeurons{1}.nIDs,mod(1:nsteps,2*trainint)==1) = HowOn; % alternating rhythmic activation of training neurons
Input(trainingNeurons{2}.nIDs,mod(1:nsteps,2*trainint)==trainint+1) = HowOn; % alternating rhythmic activation of training neurons
AltInput = Input;
trainingNeuronsAlt = trainingNeurons; clear trainingNeurons
%%

% for each seed
for seedi = 1:50
    rng(seedi);
    p.wmax = wmax; 
    p.gamma = gamma; 
    p.m = m; 
    w0 = 2*rand(p.n)*Wmax/p.n; 
    % subsong stage
    eta = p.eta;
    p.eta = 0; 
    nsteps = p.nsteps; 
    nstepsSubsong = 1000; 
    p.nsteps = nstepsSubsong; 
    w = w0; 
    p.w = w; 
    trainingNeurons{1}.nIDs = 1:k;
    trainingNeurons{2}.nIDs = 1:k;
    isOnset = rand(1,nstepsSubsong)>.9; 
    trainingNeurons{1}.tind = find(isOnset);
    trainingNeurons{2}.tind = find(isOnset); 
    trainingNeurons{1}.candLat = 1:trainint; 
    trainingNeurons{2}.candLat = 1:trainint; %1:trainint; 
    trainingNeurons{1}.thres = 12; 
    trainingNeurons{2}.thres = 12; 
    Input =-HowClamped*ones(k, nstepsSubsong); % clamp training neurons (effectively giving them higher threshold)
    Input(trainingNeurons{1}.nIDs,isOnset) = HowOn; % alternating rhythmic activation of training neurons
    Input(trainingNeurons{2}.nIDs,isOnset) = HowOn; % alternating rhythmic activation of training neurons
    bdyn = double(rand(n,nstepsSubsong)>=(1-pn)); % Random activation
    bdyn(1:k,:) = Input; 
    p.input = bdyn;
    [w xdyn] = HVCIter(p);
    LatencySubsong{seedi} = findLatency(xdyn, trainingNeurons); 
    % recovering original params
    p.eta = eta; 
    p.nsteps = nsteps; 
    
    % set(gca, 'color', 'none')
    trainingNeurons = trainingNeuronsPsyl;
    for i = 1:Niter(1)
        % Construct input
        bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
        bdyn(1:k,:) = PsylInput; 
        p.w = w; 
        p.input = bdyn;
        % One 'bout' of learning
        [w xdyn] = HVCIter(p);
    end
    LatencyEarlyPsyl{seedi} = findLatency(xdyn, trainingNeurons);
    for i = (Niter(1)+1):Niter(2)
        % Construct input
        bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
        bdyn(1:k,:) = PsylInput; 
        p.w = w; 
        p.input = bdyn;
        % One 'bout' of learning
        [w xdyn] = HVCIter(p);
    end
    
    p.wmax = wmaxSplit;
    p.m = Wmax/p.wmax;
    trainingNeurons = trainingNeuronsAlt; 
    niter = Niter(3); 
    for i = 1:niter
        % Construct input
        bdyn = double(rand(n,nsteps)>=(1-pn)); % Random activation
        bdyn(1:k,:) = AltInput; 
        p.w = w; 
        p.input = bdyn;
        p.gamma = gammas(i); 
        [w xdyn] = HVCIter(p);
    end
    LatencyLatePsyl{seedi} = findLatency(xdyn, trainingNeurons);
    seedi
end

%%
%save('C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\sigLatDistOverDev1'); 
%% compiling and plotting
%load 'C:\Users\emackev\Documents\MATLAB\code\misc_elm\HVCmodel\sigLatDistOverDev1'; 
Cases = {LatencySubsong LatencyEarlyPsyl LatencyLatePsyl}; 
figure(2);
Colors = .8*ones(3,3); 
FS = 7; % labels 
FS_axes = 5; % axis labels

for i = 1:length(Cases)
    thisCase = Cases{i};
    CompiledLats{i} = [];
    for seedi = 1:length(thisCase); 
        CompiledLats{i} = [CompiledLats{i} thisCase{seedi}{1}.mode(thisCase{seedi}{1}.FireDur)]; 
    end    
    subplot(3,1,i)
    [n,x] = hist(CompiledLats{i}, 1:10);
    bar(x*10, n/sum(n), 'facecolor', Colors(i,:))
    ylabel('Fraction','fontsize',FS)
    xlim([0 110])
    ylim([0 1.1*max(n(:))/sum(n)])
    box off
    set(gca,'color','none','tickdir','out','ticklength',[0.025 0.025])
    set(gca,'fontsize',FS_axes)
end
xlabel('Latency (ms)','fontsize',FS)


figw = 3; 
figh = 4; 
set(gcf, 'color', [1 1 1],'papersize', [figw figh], 'paperposition', [0 0 figw figh])

