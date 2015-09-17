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
