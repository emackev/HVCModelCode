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