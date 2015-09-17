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