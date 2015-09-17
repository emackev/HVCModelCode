function netForMovie_boutOnset(w, xdyn, trainingNeurons, PlottingParams)
% Makes network diagram and raster plots, called by RunHVC_boutOnset_net
% w: weight matrix
% xdyn: activity of network
% m: duration of one syllable, in timesteps
% trainingNeurons: cell array of structures containing neuron and time indices for each syllable type
% PlottingParams: sets linewidth, etc.  See RunHVC_split
%
% Emily Mackevicius 12/10/2014, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo.

%% plotting parameters

msize = PlottingParams.msize;
linewidth = PlottingParams.linewidth;
Syl1Color = PlottingParams.Syl1Color;
Syl2Color = PlottingParams.Syl2Color;
numFontSize = PlottingParams.numFontSize;
labelFontSize = PlottingParams.labelFontSize;
nplots = PlottingParams.totalPanels;
ploti = PlottingParams.thisPanel;

%% network diagram
%subplot('position', [ploti/nplots-.9/nplots, .7, .9/nplots, .2])
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
        if Latency{1}.FireDur(ni)|Latency{2}.FireDur(ni) % if it fired during either syll
            if (Latency{1}.FireDur(ni)&Latency{2}.FireDur(ni)) % if it fired during both sylls
                if (Latency{1}.mode(ni)==Latency{2}.mode(ni)) % if fired during both sylls at same phase
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
indDoubled

% keep track of which neurons participated in each syllable
FireDur1 = Latency{1}.FireDur(indkeep);
FireDur2= Latency{2}.FireDur(indkeep);

% classify neurons as specific or shared
Specific1 = FireDur1&~FireDur2;
Specific2 = FireDur2&~FireDur1;
Shared = (FireDur1&FireDur2);
indshared = find(Shared);

% calculate the incoming weights from specific neurons of each type, to
% determine sorting in y axis and color
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
% specific neurons
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

% jitter a little in x and y, so it doesn't look like a grid
jitter = .1;
% fixed seed for plotting jitter, so it doesn't interact with seed for
% running network (Seed0 = randn(1,500);)
Seed0 = [0.540299535770554,-0.332056230073109,-0.0201668962105324,1.44210050364547,-0.911471702198509,-1.27439761124562,-1.09567039916372,-0.875006453375407,1.21191939152947,0.161215819605601,0.341542595391207,-1.42947649746396,-0.210954749505668,0.341419932125287,-1.13449901886869,-0.543509331642651,-0.522115958273638,-0.0584410886801212,-0.543640928874169,0.0973942065791447,-0.729610108841670,1.53520342619730,0.181581973108703,1.42827348443680,-0.285726036162954,0.273516012508408,-1.69040905299143,-1.33896935106549,-0.755951610132419,-0.154687861692600,1.82155908115723,2.24666873905678,-0.390619025185312,0.353046852237032,-1.45891613378466,-0.109216337132524,0.876292361810659,0.155309199109896,0.229127509462572,0.857959355891777,0.311011162608834,0.107165373036265,-0.398370838058097,-0.335297040222679,0.763917702920587,1.37669204851957,-0.303884620737791,-1.50504920119019,-0.902287452911034,0.433041721281540,0.626683210517126,-0.588734983478849,-0.348207600251210,-1.01905779738314,-1.02574417480869,0.439644633444287,-0.687418468810675,-0.852618319372492,-0.846676628641204,-0.242089886164019,-0.663307090338229,-1.56674761201269,-0.0867391275060419,0.750013011258006,-1.18680682705958,0.688798255427457,0.0775804338445661,-0.408079631935128,0.666750165824117,-1.01989704905459,-0.554075735192001,-0.983310298217214,-1.77360137877207,2.53318402940323,-1.62420230399420,-0.627340386406118,-0.323222375230367,0.433099738876643,0.139204397059745,0.335297056048127,-0.257897877959559,0.498109858744791,-1.08719775857477,1.36658339036863,-1.43967307032643,0.0421470614705084,-0.677660631602860,0.0114260038561492,-0.403348428519344,0.362373763476114,-0.407115181354767,-0.926126210345786,-0.319642713927901,-1.21585416468898,-0.176003921077582,0.446128785868026,0.204469545712231,1.01026441842354,-0.243924580950639,0.471241754777045,-0.903288267236973,0.195470500004796,0.00834144483680826,-0.0567571234351842,1.14696737252755,0.481003481220535,-0.805286295816363,-1.18623215038016,-0.247198166473406,1.83929390413444,-0.0860636987504759,-1.68823071547647,0.172001147764886,1.00457559822602,0.385563747090408,1.25243652448978,-1.02295327572707,-0.813739610998985,0.552633818443304,0.168912686036287,-1.31320984077740,0.277905376984492,-0.451763100105724,1.20538229338107,0.662772428866800,-0.522370999607717,-0.363012663190467,0.355939634713330,0.206095764557025,1.56430849742699,-0.819945015599072,-0.400644621386732,0.763557158216235,-0.436189920446948,0.897736495282718,0.204495279036730,0.380321748129585,0.624645133502545,-0.100511685304420,-0.0194325044098983,1.35014348006256,-1.86456838508644,1.90559689462008,-0.114790723149111,0.126773858753657,0.545151066630993,0.661033526605043,-0.465981711452350,0.219922325382890,-0.682147520157599,-0.556805171434307,1.44111838519628,-1.25314394220787,1.73244997798176,0.151726435526522,0.407642293391344,0.820094656819088,0.151788476454889,1.14441605309776,-2.10877100039882,0.358135970067176,1.29097180157948,1.35278222286319,-0.324022750707511,1.23058622976690,-0.863424712768825,0.265097439039370,-0.256909141224611,-0.110024138278190,1.80893773641693,2.14217736400381,-0.240947268787982,-1.35173951030012,-0.818310589497384,-0.000182873609234605,-0.915857412669492,0.200859716587056,-0.293781039332376,0.648170928748985,0.0598608014404722,-0.190671912117149,2.15525801081395,-0.465917719776437,-0.121007213393658,1.52078392966238,0.928821617879422,0.357617432742821,0.228010140421689,-0.433566538110714,1.06694735296218,2.08829799214919,2.10985073820541,0.449317426324454,-0.477830642980741,-0.481983782950272,1.42970315115473,-0.178846176316974,-2.01629116266361,1.61529596998144,-0.166376100823125,1.11671490689288,1.07841012129720,-0.788274567885749,0.502292527818144,-0.321891063164824,-2.63010343559081,0.901932348613917,0.686512280198036,0.442723835746534,0.0200308821006000,-2.90646233783909,-1.50422780363728,0.160839584711452,0.461999246450974,0.254499165705607,2.02855690298655,0.975130320767919,-0.916428322939444,0.00778142012256523,-0.615203590493828,0.000122951205710950,1.13940005795281,0.243761794645890,1.24879225735953,1.42340358942916,-0.297468292862766,-0.283014798768819,0.638452404249784,0.358608358433499,-1.01444124542650,-0.545072703723796,-0.113731686968692,0.151027462113149,-1.86417992979101,-1.06473051913728,-0.853553207573177,-0.670951762961319,1.19458738584239,0.703817678016000,-0.543405241583628,-0.0575497573483473,-1.05402249602519,-0.601112719969909,-2.02167264305556,-0.464321152341921,1.73026559177343,-0.325033837429034,-0.352736411848431,-0.441212941239450,-0.105321756921205,1.29730735857052,-0.918287557738474,0.265748847598701,1.39125296758055,0.626770882775027,0.721784231460609,1.26821454496213,-0.731827229629512,0.790071873121008,-0.282432864875583,0.326652583162299,-0.363959014394143,-0.425431755856651,-1.08744480980222,-0.305381635156313,0.700757004500792,0.0490022538433330,-0.148854711439713,-1.03719679015919,-0.650837850981218,-0.599335346783658,0.535495404904010,0.260870591701303,0.0648557059925188,0.714941148645168,0.436751435691710,-1.17656951756074,0.532404016150179,-1.31042515739452,1.11260785954061,-1.00357661132050,3.18250779792889,-0.313054780742207,-0.892428687051070,0.261844895938777,-1.68598364773485,0.765635790345714,-0.664414397424032,0.895179367048265,-1.18144089872169,2.27295983567360,-0.625942500505870,-1.66324385750290,-0.672503941165623,0.877948973351547,0.567177560190190,-1.11743888366427,-0.606100876299470,-1.33598447724558,0.476224140756840,-1.19050825669453,3.26008598551012,-0.439103524847949,1.04552173622069,-0.516185377982434,-0.250202274104259,-0.0488071799936675,-0.348109367830369,-0.0120298565337762,-0.580799810702287,-0.886013436051083,0.674506620238403,-0.192611520917242,-0.981623416169810,0.816893807999206,0.00992105886174154,-0.206358323594425,0.948588631187381,-0.633137077523053,-0.208827034311530,0.0841939472841663,-0.315748270739051,0.215072050517976,1.40993922026095,1.60628930464741,1.08359143157675,-2.07959743578104,-2.23739169463690,1.49549653396230,0.572016020509796,0.0724187401439742,0.385310307138860,-1.15196646706981,-0.412444731556163,1.71117193790048,0.908866009268062,1.07735516263365,-0.906842514579264,0.658909847635298,1.27943333906709,1.05277710148249,-1.17486472195390,-0.0201403450153034,-2.44396546658145,-1.36669332949444,-1.21997642991647,-0.180384047777301,-1.42306257031966,-0.198553003553714,0.557549529845688,0.232705369184502,-0.863874906996364,-1.10441676740969,2.04215394134235,1.49208654351851,0.121554096233143,1.36083417246153,0.477772152336797,-0.763803778961587,0.204881779904442,-0.0751705201187797,0.0699001949080500,1.71246906658692,-1.25951237485211,-0.680076620143383,-0.567187534451471,0.735780752299972,0.234578810420224,1.23734421144904,-1.72913257117646,-1.40798264666776,1.63627027558892,0.852231814596294,0.333275142089968,1.18230295657475,-0.874666699691802,-0.813003479190137,-0.0317947080737779,-1.04374129490302,-1.57744156449778,0.121169744969796,-0.398445890271449,-0.118193943163204,0.429172275201091,0.536217422938857,-1.32179973239135,0.112082148856167,-1.42423816387462,1.25681933124298,0.145667848882927,-1.06255240797743,0.388985487433103,0.718754617600179,0.944193697498414,-0.505882002430843,2.41354615961092,0.0229118554189328,2.00570434487275,1.81262390284680,0.230712706802604,1.29248822018900,1.61512847428747,0.278668391679360,-0.663995267373574,-1.98867998485490,0.311485399660849,0.196831258021168,1.70402201817421,-0.226994806566507,0.239042410528162,0.415011446944945,-1.46883674420312,0.958738163001106,-2.08810552254059,-0.153692557360050,0.489669226253315,-1.22193707998467,-1.36713973643737,0.411281209612635,1.60916336125849,-0.255722787240697,-0.438399180864006,0.993290444974091,-0.146653069197029,0.378682876641687,-0.539288575143751,0.281252335653642,-0.888042200233206,0.918043721185193,0.210205395484029,0.148029593376701,0.993834934856358,-0.439151394032935,-0.349121991629529,0.542280756177789,-1.52602393404597,1.79080935300833,-1.34876705718621,0.936232038255967,-0.648450144859174,1.18545384580265,1.26483492582690,-1.32624005631073,-0.193837338107922,0.0309366260497680,0.00525415636166632,0.124395867402773,2.31339704589921,-0.209859809433100,-0.315845980384726,0.198533517186943,0.365965799503288,0.598880034580010,-0.283302185208352,1.73421549032115,-1.38221625070978,-1.24107144022877,-0.646622743668658,0.748727538378480,1.61613954732691,-1.21697934140629,0.525655224701704,0.658125949153759,-0.407773988385254,-0.896624020193576,1.16320305958340,0.637859722340896,-0.915856304833322,0.413749438211870,-0.761157093692039,0.473801592705535,-1.76847760314553,0.592399066157673,0.736053432637859,0.259312764870342,0.602169824591298,-0.524215153284913,-0.293729947956618,-0.837311466104356,0.235296118855826,-2.60437011531370,-0.565311154418695,1.86877111924631,-1.04394373191141,-0.512663549661102,-1.77850658783485,0.280439641567059,-2.05353218524893,-0.325291901441639,0.798893930832198,0.594137480832235,-1.48392866213082,0.911199553751913,-0.281483460316788,-0.854208545414047,-0.967307301023590,0.0630257497895013,-1.01255547443389,1.28829281347647,0.747933036693678;];
indJitter = setdiff(1:length(x), union(trainingset1, trainingset2)); % don't jitter seed neurons
x(indJitter)= x(indJitter)+jitter*Seed0(1:length(x(indJitter)));
y1(indJitter) = y1(indJitter)+jitter*Seed0((length(x(indJitter))+1):(2*length(x(indJitter))));

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
        loopback = (round(x(i))==round(max(x)))&(round(x(j))==round(min(x)));
        if ff & ~longrange%|loopback
            C = ones(1,3)-wplot(j,i)*ones(1,3);
            plot([x(i), x(j)], [y1(i),y1(j)], 'color', C, 'linewidth', linewidth)
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

% plot rectangle for syl1 seed neurons
rx = mean(x(trainingset1))-.5;
ry = min(y1(trainingset1))-.5;
rw = 1;
rh = max(y1(trainingset1)) - min(y1(trainingset1))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor',  'none', ...
    'LineStyle', '-', 'LineWidth', .5, 'EdgeColor', PlottingParams.SeedColor,...
    'curvature', [.98 .1])

% plot rectangle for syl2 seed neurons
rx = mean(x(trainingset2))-.5;
ry = min(y1(trainingset2))-.5;
rw = 1;
rh = max(y1(trainingset2)) - min(y1(trainingset2))+1;
rectangle('Position', [rx ry rw rh], ...
    'FaceColor',  'none', ...
    'LineStyle', '-', 'LineWidth', .5, 'EdgeColor', PlottingParams.SeedColor,...
    'curvature', [.98 .1])

xlim([trainingNeurons{1}.candLat(1)-1 trainingNeurons{1}.candLat(end)+.5])
ylim([-9 9])%[-max(abs(y1)+1) max(abs(y1)+1)]); %[min(y1)-.5 max(y1)+.5])
axis off;
set(gca, 'color', 'none')

% %% rasters...
% 
% Syl1Color = PlottingParams.Syl1Color;
% Syl2Color = PlottingParams.Syl2Color;
% ProtoSylColor = PlottingParams.ProtoSylColor;
% numFontSize = PlottingParams.numFontSize;
% labelFontSize = PlottingParams.labelFontSize;
% 
% bottom = .1;
% height = .55;
% scale = .005;
% 
% spacing = .75/(2*nplots);
% 
% % cmap = flipud(gray);
% % cmap = cmap(1:64,:);
% % cn = size(cmap,1);
% Red = trainingNeurons{1}.nIDs;
% Green = trainingNeurons{2}.nIDs;
% IsTrain1 = zeros(1,length(xdyn)); IsTrain1(Red) = 1;
% IsTrain2 = zeros(1,length(xdyn)); IsTrain2(Green) = 1;
% 
% % cmap(cn+1,:) = [0 0 0];
% % cmap(cn+2,:) = Syl1Color; % some red training neurons
% % cmap(cn+3,:) = Syl2Color; % some green training neurons
% % cmap(cn+4,:) = ProtoSylColor; % sometimes magenta ... never in this case.
% %
% % xdyn(Red,:) = xdyn(Red,:)*(1+1/cn);
% % xdyn(Green,:) = xdyn(Green,:)*(1+2/cn);
% 
% %%
% %collecting what I'll plot for the raster
% sylIDtoplot = 7; %(don't choose a protosyllable that's at the beginning of a bout)
% k = length(union(trainingset1, trainingset2));
% tindplot1 = trainingNeurons{1}.tind(sylIDtoplot) + trainingNeurons{1}.candLat-1; % time of example syl 1
% tindplot2 = trainingNeurons{2}.tind(sylIDtoplot) + trainingNeurons{2}.candLat-1; % time of example syl 2
% [~,indsort] = (sortrows(xdyn(:,[tindplot1 tindplot2]))); % sort by which fired first
% tmp = xdyn(flipud(indsort), [tindplot1 tindplot2]); % pull out the example data from xdyn
% IsTrain1 = IsTrain1(flipud(indsort));
% IsTrain2 = IsTrain2(flipud(indsort));
% indShared = (sum(tmp(:,1:length(tindplot1)),2)>0) & (sum(tmp(:,(length(tindplot1)+1):end),2)>0);
% indBO = (sum(tmp(:,1:length(tindplot1)),2)>0) & (sum(tmp(:,(length(tindplot1)+1):end),2)==0);
% rest = ~indShared;
% tmp = tmp([find(indShared); find(rest)],:); % everything that will be plotted in the rasters
% IsTrain1 = IsTrain1([find(indShared); find(rest)]);
% IsTrain2 = IsTrain2([find(indShared); find(rest)]);
% 
% 
% %%
% % plot Bout Onset syllable
% subplot('position', [ploti/nplots-2*spacing, bottom, length(tindplot1)*scale, height])%subplot(3,nHorPlot,(nHorPlot+PlottingParams.thisPanel*4-2)+[0 nHorPlot]+.75)
% tmp1 = tmp(:,1:length(tindplot1)); % just bout onset syllable
% %tmp1(end+1,end+1) = 1+4/cn; % to normalize cmap for plotting
% PlottingParams.axesPosition = [ploti/nplots-2*spacing, bottom, length(tindplot1)*scale, height];
% 
% tOffset = trainingNeurons{1}.candLat(1)-1;
% for j=1:size(tmp1,2) % for all the time steps
%     Idx = find(tmp1(1:end-1,j)>0); % find the indices of active neurons
%     if ~isempty(Idx)
%         for k=1:length(Idx) % for all the active neurons
%             Color = IsTrain1(Idx(k))*PlottingParams.Syl1Color + ...
%                 IsTrain2(Idx(k))*PlottingParams.Syl2Color;
%             h = patch(10*([j-1,j,j,j-1]+tOffset),[Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],Color,'edgecolor','none');
%         end
%     end
% end
% 
% hold on; box off
% set(gca, 'fontsize', numFontSize)
% set(gca, 'color', 'none', 'xtick', [0 50 100], 'xticklabel', {'0', '50', '100'},'ydir', 'reverse', 'fontsize', numFontSize)
% set(gca, 'ydir', 'reverse','tickdir','out','ticklength',[0.015 0.015], 'color', 'none', 'fontsize', numFontSize,'tickdir','out');
% 
% if PlottingParams.thisPanel==1
%     ylabel('Neuron','fontsize', labelFontSize)
%     set(gca,'ytick',0:20:100,'fontsize', numFontSize)
% else
%     set(gca,'ytick',0:20:100,'yticklabel',{});
% end
% 
% if PlottingParams.Hor
%     if sum(indShared)>0 % if shared neurons
%         plot([-10+trainingNeurons{1}.candLat(1)*10 trainingNeurons{1}.candLat(end)*10], (sum(indShared))*ones(1,2), 'k', 'linewidth', PlottingParams.linewidth);
%     end
%     plot([-10+trainingNeurons{1}.candLat(1)*10 trainingNeurons{1}.candLat(end)*10], ...
%         (sum(indBO)+sum(indShared))*ones(1,2), 'k', 'linewidth', PlottingParams.linewidth);
% end
% if isfield(PlottingParams,'boutOnsetElement')
%     patch([-10+trainingNeurons{1}.candLat(1)*10 -10 -10 -10+trainingNeurons{1}.candLat(1)*10],[-4 -4 -2 -2],Syl1Color);
%     patch([-10+trainingNeurons{2}.candLat(1)*10 trainingNeurons{2}.candLat(end)*10 trainingNeurons{2}.candLat(end)*10 -10+trainingNeurons{2}.candLat(1)*10],[-4 -4 -2 -2],Syl2Color);
%     if PlottingParams.thisPanel>1 % TO
%         text((-10+trainingNeurons{1}.candLat(1)*10-10)/2,-7,'\epsilon','fontsize',7)
%         text(((-10+trainingNeurons{2}.candLat(1)*10)+(trainingNeurons{2}.candLat(end)*10))/2,-7,'\alpha','fontsize',7)
%     end
% else
%     patch([-10+trainingNeurons{1}.candLat(1)*10 trainingNeurons{1}.candLat(end)*10 trainingNeurons{1}.candLat(end)*10 -10+trainingNeurons{1}.candLat(1)*10],[-4 -4 -2 -2],Syl1Color);
%     if PlottingParams.thisPanel>2
%         text(((-10+trainingNeurons{1}.candLat(1)*10)+(trainingNeurons{1}.candLat(end)*10))/2,-7,'\beta','fontsize',7)
%     end
% end
% ylim([-5 ntot])
% xlim([-10+trainingNeurons{1}.candLat(1)*10 trainingNeurons{1}.candLat(end)*10+10])
% %%
% 
% % Plot protosyllable.
% subplot('position', [ploti/nplots-spacing, bottom, length(tindplot2)*scale, height])%subplot(3,nHorPlot,(nHorPlot+PlottingParams.thisPanel*4)+[0 nHorPlot])
% tmp1 = tmp(:,(length(tindplot1)+1):end); % just for protosyllable
% %tmp1(end+1,end+1) = 1+4/cn; % to normalize cmap for plotting
% PlottingParams.axesPosition = [ploti/nplots-spacing, bottom, length(tindplot2)*scale, height];
% 
% tOffset = trainingNeurons{2}.candLat(1)-1;
% for j=1:size(tmp1,2) % for all the time steps
%     Idx = find(tmp1(1:end-1,j)>0); % find the indices of active neurons
%     if ~isempty(Idx)
%         for k=1:length(Idx) % for all the active neurons
%             Color = IsTrain1(Idx(k))*PlottingParams.Syl1Color + ...
%                 IsTrain2(Idx(k))*PlottingParams.Syl2Color;
%             h = patch(10*([j-1,j,j,j-1]+tOffset),[Idx(k)-1,Idx(k)-1,Idx(k),Idx(k)],Color,'edgecolor','none');
%         end
%     end
% end
% 
% % PlottingParams.tOffset = trainingNeurons{2}.candLat(1)-1;
% % plotRaster(tmp1,PlottingParams,0,2)
% %imagesc(tmp1, 'xdata', trainingNeurons{2}.candLat*10); colormap(cmap); axis tight
% hold on; box off
% set(gca, 'fontsize', numFontSize)
% set(gca, 'color', 'none', 'xtick', [0 50 100], 'xticklabel', {'0', '50', '100'}, 'ydir', 'reverse', 'fontsize', numFontSize)
% set(gca, 'ydir', 'reverse','tickdir','out','ticklength',[0.015 0.015], 'color', 'none', 'fontsize', numFontSize,'tickdir','out');
% set(gca, 'ytick',0:20:100,'yticklabel',{})
% if PlottingParams.Hor
%     if sum(indShared)>0 % if shared neurons
%         plot([-10+trainingNeurons{2}.candLat(1)*10 trainingNeurons{2}.candLat(end)*10], (sum(indShared))*ones(1,2), 'k', 'linewidth', PlottingParams.linewidth);
%     end
%     plot([-10+trainingNeurons{2}.candLat(1)*10 trainingNeurons{2}.candLat(end)*10], ...
%         (sum(indBO)+sum(indShared))*ones(1,2), 'k', 'linewidth', PlottingParams.linewidth);
% end
% 
% patch([-10+trainingNeurons{2}.candLat(1)*10 trainingNeurons{2}.candLat(end)*10 trainingNeurons{2}.candLat(end)*10 -10+trainingNeurons{2}.candLat(1)*10],[-4 -4 -2 -2],Syl2Color);
% if isfield(PlottingParams,'boutOnsetElement')
%     if PlottingParams.thisPanel>1
%         text(((-10+trainingNeurons{2}.candLat(1)*10)+(trainingNeurons{2}.candLat(end)*10))/2,-7,'\alpha','fontsize',7)
%     end
% else
%     if PlottingParams.thisPanel>2
%         text(((-10+trainingNeurons{2}.candLat(1)*10)+(trainingNeurons{2}.candLat(end)*10))/2,-7,'\alpha','fontsize',7)
%     end
% end
% ylim([-5 ntot])
% xlim([-10+trainingNeurons{2}.candLat(1)*10 trainingNeurons{2}.candLat(end)*10+10])