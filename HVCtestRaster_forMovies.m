function HVCtestRaster_forMovies(xdyn,Input,w, plottingParams)

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
cmap = jet;
cmap(1,:) = ones(1,3);
%cmap(2,:) = .5*ones(1,3); 

% h = subplot(2,2,1);cla;  
% hold on; set(gca, 'ydir', 'reverse'); xlim([0 10]); ylim([0 size(w,1)])
% lRaster = 10; 
% aAligned = zeros(size(xdyn,1),lRaster+1);
% for i = 1:length(indA)
%     if (indA(i)+lRaster)<=size(xdyn,2)
%         aAligned = aAligned + xdyn(:,indA(i):(indA(i)+lRaster)); 
%     end
% end
% [~,sortInda] = sortrows(aAligned((k+1):end,:)); 
% sortInda = [1:k k+flipud(sortInda)'];
% [~,unsorta] = sort(sortInda); 
% aAligned = aAligned(sortInda,:); 
% aAligned = aAligned>thres;aAlignedBin=aAligned;
% aAligned = bsxfun(@times, aAligned, (1:size(aAligned,1))');
% aAligned(1) = size(aAligned,1); % to scale colormap
% imagesc(aAligned, 'xdata', (1:(lRaster+1))); colormap(cmap); 
% 
% g = subplot(2,2,3); cla; 
% grid on; set(gca, 'ydir', 'reverse'); xlim([0 10]); ylim([0 size(w,1)])
% lRaster = 10; 
% bAligned = zeros(size(xdyn,1),lRaster+1);
% for i = 1:length(indB)
%     if (indB(i)+lRaster)<=size(xdyn,2)
%         bAligned = bAligned + xdyn(:,indB(i):(indB(i)+lRaster)); 
%     end
% end
% [~,sortIndb] = sortrows(bAligned((k+1):end,:)); 
% sortIndb = [1:k k+flipud(sortIndb)'];
% [~,unsortb] = sort(sortIndb); 
% bAligned = bAligned(sortIndb,:); 
% bAligned = bAligned>thres;
% colorInd = 1:size(aAligned,1); 
% colorInd(sum(aAlignedBin,2)<=thres) = 2; % color specific neurons gray
% colorInd = colorInd(unsorta); % in original sorting 
% colorInd =  colorInd(sortIndb);
% bAligned = bsxfun(@times, bAligned, colorInd');
% bAligned(1) = size(aAligned,1); % to scale colormap
% imagesc(bAligned, 'xdata', (1:lRaster+1)-1); colormap(cmap);  
% linkaxes([h g], 'x'); 

subplot(3,plottingParams.totalPanels,plottingParams.thisPanel)
cla; 
k = size(Input,1);
tind = 1:80; 
bOnOffset = diff(find(sum(Input,1)));
trainint = bOnOffset(2);
bOnOffset = bOnOffset(1);
n = size(xdyn,1); 
set(gca, 'color', 'none')
xplot = xdyn(:, tind);
xplot = xplot+repmat(.05*(sum(Input(:,1:size(xplot,2))>0,1)>0), size(xplot,1),1);
% sortFrom = find(Input(end,:)>0); sortFrom = sortFrom(1); 
% [~,sortInd] = sortrows(xdyn((k+1):n,sortFrom:end)); 
sortInd = [1:(k-1) flipud(sortbyCorr(w(k:end,k:end))+k-1)]; 
xplot = xplot((sortInd),tind);
% [~,sortInd2] = sort(sum(xplot,2),'descend'); 
% xplot = xplot(sortInd2,:);
imagesc([(1:size(xdyn,1))'/size(xdyn,1) xplot]); colormap(cmap)%(flipud(gray))
xlabel('Time (10 ms)', 'fontsize', labelFontSize)
ylabel('Neuron', 'fontsize', labelFontSize)
set(gca, 'fontsize', numFontSize)


subplot(3,plottingParams.totalPanels,plottingParams.totalPanels+plottingParams.thisPanel)
imagesc(w(sortInd,sortInd))
ylabel('Neuron', 'fontsize', labelFontSize)
xlabel('Neuron', 'fontsize', labelFontSize)
set(gca, 'fontsize', numFontSize)
%%
subplot(3,plottingParams.totalPanels,2*plottingParams.totalPanels+plottingParams.thisPanel); cla; hold on
%set(gca, 'color', zeros(1,3))
w1 = w(sortInd,sortInd);
n = size(w,1); 
S = w1+w1'+(w1*w1')/max(max(w1*w1')); 
S = S/max(S(:));
%S = bsxfun(@rdivide, S, sum(S.*eye(n)));
S = S.*~eye(n)+eye(n); 
S = conv2(S, gausswin(5)*gausswin(5)', 'same'); 
D = (1-S); 
D = D-min(D(:)); 
D = D.*~eye(size(w,1)); 
D = D/max(D(:));
[Y,e] = cmdscale(D);
%Y = mdscale(D,2);%, 'Criterion', 'strain');
Colors = jet(size(w,1)); 
% for i = 1:size(w,1)
%     for j = i:size(w,2)
%         ind = [i j];
%         if w1(i,j)>prctile(w1(:), 0)
%             plot(Y(ind,1), Y(ind,2),'k', 'color', w1(j,i)/max(w1(:))*ones(1,3))%, 'linewidth', lwidth, 'Markersize', msize);shg
%         end
%     end
% end
axis off
scatter(Y(:,1), Y(:,2), '.', 'cdata', Colors, 'sizedata', 1000*ones(1,n))