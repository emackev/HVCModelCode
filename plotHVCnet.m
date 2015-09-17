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
