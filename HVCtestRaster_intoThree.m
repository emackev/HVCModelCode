function HVCtestRaster_intoThree(xdyn,Input,w, plottingParams)
%%
if nargin < 4
    plottingParams.totalPanels = 1;
    plottingParams.thisPanel = 1; 
    plottingParams.numFontSize = 5; 
    plottingParams.labelFontSize = 6; 
end;

numFontSize = plottingParams.numFontSize;
labelFontSize = plottingParams.labelFontSize;


k = size(Input,1);
thres = 1; 

indA = find(Input(1,:));
indB = find(Input(end,:));
cmap = flipud(hot);
cmap(1,:) = ones(1,3);

subplot(1,plottingParams.totalPanels,plottingParams.thisPanel)
cla; 
k = size(Input,1);
tind = 1:60; 
bOnOffset = diff(find(sum(Input,1)));
trainint = bOnOffset(2);
bOnOffset = bOnOffset(1);
n = size(xdyn,1); 
set(gca, 'color', 'none')
xplot = xdyn(:, tind);
%xplot = xplot+repmat(.2*mod(cumsum(sum(Input(:,1:size(xplot,2))>0,1)>0)+2,3), size(xplot,1),1);
xplot(xplot>1)=1; 
if issame(plottingParams.sortby, 'activity')
    sortFrom = find(Input(1,:)>0); sortFrom = sortFrom(1); 
    [~,sortInd] = sortrows(xdyn((k+1):n,sortFrom:end)); 
    sortInd = [(1:k)'; k+flipud(sortInd)];
else % sort by weight matrix
    sortInd = [1:(k-1) flipud(sortbyCorr(w(k:end,k:end))+k-1)]; 
end
xplot = xplot(sortInd,:);
%
if ~issame(Input(1,:), Input(end,:)) % if in splitting phase
    indSharedThree = sum(xplot((k+1):n,:)==1,2)>4;
    indSharedTwo = (sum(xplot((k+1):n,:)==1,2)>2)&(sum(xplot((k+1):n,:)==1,2)<5); 
    rest = (~indSharedThree)&(~indSharedTwo);
    sortInd = ([(1:k)'; find(indSharedThree)+k; find(indSharedTwo)+k; find(rest)+k]);
    xplot = xplot(sortInd,:); 
    isProto = zeros(1,size(xplot,1)); 
    isSyl1 = zeros(1,size(xplot,1)); isSyl1(1:3) = 1; 
    isSyl2 = zeros(1,size(xplot,1)); isSyl2(4:6) = 1; 
    isSyl3 = zeros(1,size(xplot,1)); isSyl3(7:9) = 1; 
else
    isProto = zeros(1,size(xplot,1)); isProto(1:k) = 1; 
    isSyl1 = zeros(1,size(xplot,1));
    isSyl2 = zeros(1,size(xplot,1));
    isSyl3 = zeros(1,size(xplot,1));
end
xplot = xplot(:,tind);
tOffset = 0; 
psylColor = [1 0 1]; 
for j=1:size(xplot,2) % for all the time steps
    Idx = find(xplot(1:end-1,j)>0); % find the indices of active neurons    
    if ~isempty(Idx)
        for k=1:length(Idx) % for all the active neurons
            Color = isProto(Idx(k))*psylColor + ...
                isSyl1(Idx(k))*[0 0 1] + ...
                isSyl2(Idx(k))*[0 1 1] + ...
                isSyl3(Idx(k))*[1 0 0];
            h = patch(10*([j-1,j,j,j-1]+tOffset),[Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],Color,'edgecolor','none');
        end  
    end
end

set(gca, 'ydir', 'reverse')
k = size(Input,1);
%imagesc(xplot, 'xdata',tind*10); colormap(cmap)%(flipud(gray))
xlabel('Time (ms)', 'fontsize', labelFontSize)

if plottingParams.thisPanel==1
    ylabel('Neuron', 'fontsize', labelFontSize)
end
hold on
plot([0 size(xplot,2)*10], (k)*ones(1,2), 'k', 'linewidth',  plottingParams.linewidth); 

if ~issame(Input(1,:), Input(end,:)) % if in splitting phase
    plot([0 size(xplot,2)*10], (k+sum(indSharedThree)+sum(indSharedTwo))*ones(1,2), 'k', 'linewidth', plottingParams.linewidth); 
    sylStarts = 0:100:500; 
    for i = 1:6
        Color = (mod(i,3)==0)*[1 0 0] + (mod(i,3)==1)*[0 0 1] + (mod(i,3)==2)*[0 1 1];
        patch([0 80 80 0]+sylStarts(i),[-4 -4 -2 -2],Color);
        
        %% TO
        switch mod(i,3)
            case 1
                text(25+sylStarts(i),-6,'A','fontsize',7)
            case 2
                text(25+sylStarts(i),-6,'B','fontsize',7)
            case 0
                text(25+sylStarts(i),-6,'C','fontsize',7)
        end
    end
else
    sylStarts = 0:100:500; 
    for i = 1:6
        patch([0 80 80 0]+sylStarts(i),[-4 -4 -2 -2],[1 0 1]);
    end
end
set(gca, 'fontsize', numFontSize)
axis tight
ylim([-5 size(xdyn,1)]); 
set(gca,'tickdir','out','ticklength',[0.015 0.015], 'color', 'none', 'fontsize', numFontSize);

if plottingParams.thisPanel==1
    set(gca,'ytick',0:20:100)
else
    set(gca,'ytick',0:20:100,'yticklabel',{})
end

if plottingParams.thisPanel == 4
    tStart = [1 10 20]; 
    for i = 1:2
        spec = find(sum(xplot(:,(tStart(i)+1):(tStart(i)+10)),2)); 
        %spec = spec(spec>9); 
        spec = spec(end); 
        plot([0 size(xplot,2)*10], [spec spec], 'k',  'linewidth', plottingParams.linewidth);
    end
end



% subplot(2,plottingParams.totalPanels,plottingParams.totalPanels+plottingParams.thisPanel)
% imagesc(w(sortInd,sortInd))
% ylabel('Neuron', 'fontsize', labelFontSize)
% xlabel('Neuron', 'fontsize', labelFontSize)
% set(gca, 'fontsize', numFontSize)