% VERSION HISTORY
%
% 231115 : Idse - change to deal with self-stitched files
%

clear; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);
% 
% % define channels and conditions
% %--------------------------------------
% manMeta = struct();
% manMeta.nucChannel = 0;
% manMeta.channelLabel = {'DAPI','BCAT','BRA','pFGFR1'};
% manMeta.conditions = {'B50 FGFRi', 'B50', 'mTeSR'};
% manMeta.nChannels = numel(manMeta.channelLabel);
% 
% % proper preprocessing should make it unnecessary to define the resolution
% manMeta.zres = 2;
% manMeta.xres = 0.353;
% manMeta.yres = 0.353;
% 
% % define micropatterned 'MP' or 'Disordered'
% imageType = 'MP'; 
% % radii of micropatterns in micron for each condition (ignored if not 'MP')
% radii = [700 700 700]/2; 
% % margin in microns to keep around the nominal radius (for masking cells)
% Rmarg=55; 
% 
% % different cases for input files
% %--------------------------------------
% 
% % self-stitched - channels are split
% if ~isempty(dir(fullfile(dataDir,'stitched*.tif')))
%     filelistAll = {};
%     filenameFormat = "stitched_p%.4d_w%.4d_t%.4d.tif";
%     xsectfname_prefix_format = 'stitched_p%.4d';
%     
% % Dragonfly stitched - channels are not split
% elseif ~isempty(dir(fullfile(dataDir,'*FusionStitcher*.ims')))
%     filelistAll = dir(fullfile(dataDir,'*FusionStitcher*.ims'));
%     % FINISH
% end
% 
% manMeta.nWells = numel(manMeta.conditions);
% manMeta.posPerCondition = [4 5 4];
% manMeta.nPositions = sum(manMeta.posPerCondition);

% reload previous metadata or create new
%--------------------------------------
metafile = fullfile(dataDir,'meta.mat');
if exist(metafile)
    load(metafile,'meta');
    disp(['loading meta file ' metafile]);
else
    meta = Metadata(dataDir, manMeta);
    save(metafile, 'meta');
    disp(['created meta file ' metafile]);
end

% reload previous results
%--------------------------------------
posfile = fullfile(dataDir,'positions.mat');
if exist(posfile)
    load(posfile, 'positions');
    disp(['loading positions file ' posfile]);
end

statsfile = fullfile(dataDir,'stats.mat');
if exist(statsfile)
    load(statsfile, 'stats');
    disp(['loading stats file ' statsfile]);
end

countfile = fullfile(dataDir,'counts.mat');
if exist(countfile)
    load(countfile, 'counts');
    disp(['loading count file ' countfile]);
end


% ratio between z and xy resolution
scalefactor = meta.zres / meta.xres;

%% cleanup stack for specific positions

zcorrection = true;
zcorrfactor = [1.02 1.02 1 1];
mediansize = [1 1];

tic
%parpool(3)

for pi = [2 8]%[2 5 8]
    
        fname = char(sprintf(positions(pi).filenameFormat, pi-1));
        cleanupStack(dataDir, fname, zcorrection, zcorrfactor, mediansize);
end
toc


%%
%--------------------------------------------------------------------------
% INITIAL QC
%--------------------------------------------------------------------------

%% overlay quantification for QC on XY MIP

piSerious = 2;%[1,5,9];
ci = 3; % define nucChannel you want to check
checkOPT = 'nuc'; % choose from 'nuc' 'cyto' 'NCratio'

for pi = piSerious
    condi = find(pi >= meta.conditionStartPos,1,'last');

    tol = 0.01;
    img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);

    % saturate ALL PIXEL according to tol
    imgGREY = imadjust(img,stretchlim(img,[tol,1-tol]),[]);
    imgBLACK = imbinarize(img, 1);

    nucLevel = positions(pi).cellData.nucLevel;
    cytLevel = positions(pi).cellData.cytLevel;
    background = positions(pi).cellData.background;
    XY = positions(pi).cellData.XY;

    cNuc = nucLevel(:,ci)-background(ci);
    cCyt = cytLevel(:,ci)-background(ci);
    NC = cNuc./cCyt;

    if strcmp(checkOPT,'nuc')
        z = cNuc;
    elseif strcmp(checkOPT,'cyto')
        z = cCyt;
    elseif strcmp(checkOPT,'NCratio')
        z = NC;
    end

    % saturate ALL CELLS according to tol
    n = length(z);
    zs = sort(z);
    zmin = zs(max(ceil(n*tol),1));
    zmax = zs(floor(n*(1-tol)));
    z(z < zmin) = zmin; z(z > zmax) = zmax;
    z = (z-min(z))/(max(z)-min(z));

    figure()
    imshow(imgGREY)
    hold on
    s = scatter(XY(:,1), XY(:,2), 50, z*[1, 1, 1], "filled");
    hold off
    %s.MarkerEdgeColor = 'm';
    title(['Intensity in ',meta.channelLabel{ci}])
    saveas(gcf,['Intensity_check_' meta.channelLabel{ci} '_' meta.conditions{condi} '_Pos' num2str(pi) '.jpg']);
end




%% get stats

stats = cellStats(positions, meta, positions(1).dataChannels);
% automaticly get thresholds
confidence = 0.95;
conditions = 1:numel(meta.conditions);
whichthreshold = []; %[1 1 2]
stats.getThresholds(confidence, conditions, whichthreshold);

% Make nucHistogram channel by channel
conditionIdx = 1:numel(meta.conditions); % [1,3,5] Specify conditions for Histograms
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
for channelIdx = 1:numel(manMeta.channelLabel)
    options = struct('channelIndex',channelIdx, 'cytoplasmic', false,...
        'cumulative',false, 'time', 1,...
        'conditionIdx',conditionIdx,...
        'titlestr',meta.channelLabel{channelIdx}...
        );
    figure()
    stats.plotDistributionComparison(options)
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_dist.png']));
end
%ylim([0 0.02])
%xlim([0 1000])

save(statsfile, 'stats');
%stats.exportCSV(meta); %export to CSV

%% check if scatter plots of markers in different samples in same conditions look similar

condi = 1; % define the condition you want to check
i = 2; % define two channels you want to compare
j = 3; % define two channels you want to compare

X = stats.nucLevel{condi}(:,i);
Y = stats.nucLevel{condi}(:,j);
X = log(1+X/mean(X));
Y = log(1+Y/mean(Y));
scatter(X, Y, 10, stats.sample{condi}, 'filled')
% xlim([0 1.5]);
% ylim([0 3]);

%% Create combined nucHistogram channel by channel, If mannually thresholding, use this as a referrence

nucLevelCombined = [];
XYCombined = [];
sampleIdx = [];
for i = 1:numel(manMeta.conditions)
    nucLevel = stats.nucLevel{i};
    nucLevelCombined = [nucLevelCombined; nucLevel];
end
for i=1:numel(manMeta.channelLabel)
    figure()
    histogram(nucLevelCombined(:,i));
    xlabel('Intensity');
    ylabel('Population');
    title([manMeta.channelLabel{i},' All'])
    %     xlim([0,3000]);
    %     ylim([0,3000]);
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{i} '_tot_dist.jpg']));
end


%% NORMALIZE BY DAPI FOR INTENSITY VARIATION BETWEEN WELLS
% 
% statsNorm = load(statsfile, 'stats');
% statsNorm = statsNorm.stats;
% 
% ci = 1;
% meannucI = zeros([numel(meta.conditions) 2]);
% 
% for condi = 1:numel(meta.conditions)
%         meannucI(condi, 1) = condi;
%         meannucI(condi, 2) = mean(statsNorm.nucLevel{condi}(:,ci));
% end
% % scatter(meannucI(:,1),meannucI(:,2))
% % ylim([0 max(meannucI(:,2))])
% 
% N = mean(meannucI(:,2));
% Iscale = zeros([numel(meta.conditions) 1]);
% 
% for ci = 2:4
% 
%     meanI = zeros([numel(meta.conditions) 2]);
%     for condi = 1:numel(meta.conditions)
% 
%             Iscale(condi) = N/mean(statsNorm.nucLevel{condi}(:,1));
% 
%             meanI(condi, 1) = condi;
%             meanI(condi, 2) = Iscale(condi)*mean(statsNorm.nucLevel{condi}(:,ci));
%             statsNorm.nucLevel{condi}(:,ci) = Iscale(condi)*statsNorm.nucLevel{condi}(:,ci);
%             statsNorm.cytLevel{condi}(:,ci) = Iscale(condi)*statsNorm.cytLevel{condi}(:,ci);
%     end
% 
% end
% %scatter(meanI(:,1),meanI(:,2))

%% make intensity radial profile (MP Only)

% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', method 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
options = struct('nucChannels', [2], ...
    'cytChannels', [3],...'ratioMode', 'N:C',...
    'normalize', true, 'std', true, 'FontSize', 15, 'legend',false,...
    'pointsPerBin', 300, 'interpmethod', 'linear', 'stdMethod', 'perColony',...
    'log1p',false);

% options = struct('ratioChannels', 1, 'ratioMode', 'N:C',...
%     'normalize', true, 'std', true, 'FontSize', 15, 'legend',true,...
%     'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
%     'log1p',false);

options.colors = lines(5);
options.colors = options.colors([5 2 1 3 4],:);

res = {};
% first round without normalization (to collect overall max and min of all), keep normalize = false

conditionindices = 1:3;%[13 16:18];
%conditionindices = [13:15];

for condi = conditionindices%1:numel(meta.conditions)
    clf
    
    tempoptions = options;
    tempoptions.normalize = false;
    tempoptions.conditionIdx = condi;
    tempoptions.colonyRadius = radii(condi);
    
    res{condi} = plotRadialProfiles(stats, meta, tempoptions);
    
    min(res{condi}.ratio_profile)
    
    if condi == conditionindices(1)
        maxvalsNucAll = max(res{condi}.nuc_profile);
        minvalsNucAll = min(res{condi}.nuc_profile);

        maxvalsCytAll = max(res{condi}.cyt_profile);
        minvalsCytAll = min(res{condi}.cyt_profile);

        maxvalsRatioAll = max(res{condi}.ratio_profile);
        minvalsRatioAll = min(res{condi}.ratio_profile);
    else
        maxvalsNucAll = max(max(res{condi}.nuc_profile), maxvalsNucAll);
        minvalsNucAll = min(min(res{condi}.nuc_profile), minvalsNucAll);

        maxvalsCytAll = max(max(res{condi}.cyt_profile), maxvalsCytAll);
        minvalsCytAll = min(min(res{condi}.cyt_profile), minvalsCytAll);

        maxvalsRatioAll = max(max(res{condi}.ratio_profile), maxvalsRatioAll);
        minvalsRatioAll = min(min(res{condi}.ratio_profile), minvalsRatioAll);
    end
end
close all

% second round with normalization
nuclimits = [minvalsNucAll' maxvalsNucAll'];
options.nuclimits = nuclimits;

cytlimits = [minvalsCytAll' maxvalsCytAll'];
options.cytlimits = cytlimits;

ratiolimits = [minvalsRatioAll' maxvalsRatioAll'];
options.ratiolimits = ratiolimits;

for condi = conditionindices
    lw = 5;
    options.lw = lw;
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    %options.legend = false;
    figure()
    res{condi} = plotRadialProfiles(stats, meta, options);
    ylim([0 1.1]);
    ylabel('')
    %ylim([-0.0195 -0.016]);
    xlim([0 radii(condi)]);
    fs = 36;
    axis square;
    cleanSubplot(fs,lw)
    %title(meta.conditions{condi})
    saveas(gcf, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.png'])
    %writetable(res{condi}.datatable, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.csv'])
    close;
end

%% radial intensity - conditions combined

markerChannels = [2 3 4];%options.nucChannels+1;
conditionindices = [2 1];
legendstr = {'BMP4','+FGFRi'};

for channel = markerChannels
    
    figure()
    colors = lines(numel(conditionindices));

    lw = 3;
    fs = 30;
    lfs = 24;
    %legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.nuc_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        %legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
    end
    i = 0;
    for condi = conditionindices
        i = i+1;
        
        plerr = res{condi}.nuc_profile(:,channel) + res{condi}.nuc_profile_colstd(:,channel);
        minerr = res{condi}.nuc_profile(:,channel) - res{condi}.nuc_profile_colstd(:,channel);
        good = ~isnan(plerr);
        fill([res{condi}.r(good),fliplr(res{condi}.r(good))],[plerr(good)', fliplr(minerr(good)')],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    hold off
    xlim([0 350]);
    ylim([0 1.25]);
    xlabel('edge distance ( um )')
    if options.log1p
        ylabel('log( intensity )');
    else
        ylabel('intensity');
    end
    legend(p,legendstr,'FontSize',lfs,'Location','NorthWest')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    axis square;
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.png']));
    close;
end

%% radial intensity - conditions combined for ratio

markerChannels = options.ratioChannels+1;

for channel = markerChannels
    figure()
    colors = lines(numel(conditionindices));

    lw = 2;
    fs = 26;
    legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.ratio_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
    end
    i = 0;
    for condi = conditionindices
        i = i+1;
        
        plerr = res{condi}.ratio_profile(:,channel) + res{condi}.ratio_profile_std(:,channel);
        minerr = res{condi}.ratio_profile(:,channel) - res{condi}.ratio_profile_std(:,channel);
        good = ~isnan(plerr);
        fill([res{condi}.r(good),fliplr(res{condi}.r(good))],[plerr(good)', fliplr(minerr(good)')],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    hold off
    xlim([0 350]);
    ylim([-0.1 1.5]);
    xlabel('edge distance ( um )')
    if options.log1p
        ylabel('log( intensity )');
    else
        ylabel('N:C intensity');
    end
    legend(p,legendstr,'FontSize',15,'Location','best')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    axis square;
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.png']));
end

%% make intensity radial profile, one colony (MP Only)

% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
% 4. decrease pointsPerBin when you analysis single colony
posList = 1:13;
norMethodPos = 1; % 1 - normalize individually; 2 - normalize according to all position
options = struct('nucChannels', 1:3, 'ratioMode', 'C:N',...
    'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 100, 'interpmethod', 'nearest', 'stdMethod', 'perColony');
% first round without normalization (to collect overall max and min of all), keep normalize = false
count = 0;
for posidx = posList
    
    condi = find(pi >= meta.conditionStartPos,1,'last');
    resOP = {};
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    metaOP.posPerCondition = 1;

    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);

    clf
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    close all

    count = count+1;
    if count == 1
        maxvalsNuc = max(resOP.nuc_profile);
        minvalsNuc = min(resOP.nuc_profile);

        maxvalsCyt = max(resOP.cyt_profile);
        minvalsCyt = min(resOP.cyt_profile);

        maxvalsRatio = max(resOP.ratio_profile);
        minvalsRatio = min(resOP.ratio_profile);
    else
        maxvalsNuc = max(maxvalsNuc, max(resOP.nuc_profile));
        minvalsNuc = min(minvalsNuc, min(resOP.nuc_profile));

        maxvalsCyt = max(maxvalsCyt, max(resOP.cyt_profile));
        minvalsCyt = min(minvalsCyt, min(resOP.cyt_profile));

        maxvalsRatio = max(maxvalsRatio, max(resOP.ratio_profile));
        minvalsRatio = min(minvalsRatio, min(resOP.ratio_profile));
    end
end

if norMethodPos == 2
    nuclimits = [minvalsNuc' maxvalsNuc'];
    options.nuclimits = nuclimits;
    cytlimits = [minvalsCyt' maxvalsCyt'];
    options.cytlimits = cytlimits;
    ratiolimits = [minvalsRatio' maxvalsRatio'];
    options.ratiolimits = ratiolimits;
end

options.normalize = true; % set normalization here

for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition(1));
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    metaOP.posPerCondition = 1;
    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);
    figure()
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    ylim([0 1.1]);
    xlim([0 radii(condi)]);
    axis square
    title([meta.conditions{condi} ' Pos' num2str(posidx)])
    saveas(gcf, ['radialProfileOP_' meta.channelLabel{:} '_' meta.conditions{condi} '_Pos' num2str(posidx) '.png'])
end

%% scatter plot subpopulation on XY MIP to check thresholds

condi = 1; % condition index
cpi = 3; % position within condition index
s = 30;

pi = cpi + meta.conditionStartPos(condi);
%subDir = filelist{condi}(pi).folder;

channelThresholds = [0, stats.thresholds(2), stats.thresholds(3), stats.thresholds(4)]; % use automatic thresholds
%channelThresholds = [0 500 500 400]; % manually adjust thresholds

nucLevel = positions(pi).cellData.nucLevel;
background = positions(pi).cellData.background;
positive = {};
for i = 1:meta.nChannels
    positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
end

for ci = 2%:numel(manMeta.channelLabel)
    
    % s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
    % prefix = [s{:}];
    %
    % zdir = [prefix '_zslices'];
    % img = imread(fullfile(subDir, zdir, sprintf([prefix '_MIP_w%.4d.tif'], ci-1)));
    
    img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);
    img = imadjust(img,stitchedlim(img));
    %img = mat2gray(img, [150 4000]);
    [X,Y] = meshgrid(1:size(img,2),1:size(img,1));
    if strcmp(imageType,'MP')
        R = sqrt((X - positions(pi).center(1)).^2 + (Y - positions(pi).center(2)).^2);
        mask = R > positions(pi).radiusPixel + Rmarg/meta.xres;
        img(mask) = max(img(:));
    end

    figure,
    imshow(img)
    hold on
    %scatter(positions(pi).cellData.XY(:,1),positions(pi).cellData.XY(:,2),'filled')
    XY = positions(pi).cellData.XY(positive{ci},:);
    scatter(XY(:,1),XY(:,2),s,'filled','r')
    %XY = positions(pi).cellData.XY(AP2Cp,:);
    %scatter(XY(:,1),XY(:,2),'x','g')
    XY = positions(pi).cellData.XY(positive{2},:);
    %scatter(XY(:,1),XY(:,2),'x','b')
    if strcmp(imageType,'MP')
        scatter(positions(pi).center(1), positions(pi).center(2),500,'.','g')
    end
    hold off
    title(meta.channelLabel{ci})
    saveas(gcf,['subpopulation_' meta.channelLabel{ci} '_Pos' num2str(pi) '.jpg']);
end
% save new manually adjusted thresholds
for pi = 1:numel(positions)
    positions(pi).cellData.channelThresholds = channelThresholds;
end
save(fullfile(dataDir,'positions'), 'positions');
stats.thresholds = channelThresholds;
save(statsfile, 'stats');


%% radial probability of being positive by condition (MP Only)

Pall = {};
Pstdall = {};
xall = {};

options = struct('edgeDistance',false);

for condi = 1:numel(meta.conditions)
    combos = {}; % [2 4],[2 3],[3 4]};
    stats.markerChannels = 2:4;
    [P,x,Pstd] = radialPositive(stats, condi, meta, combos, options);
    Pall{condi} = P;
    Pstdall{condi} = Pstd;
    xall{condi} = x;
    ylim([0 1]);
    saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{condi} '_', meta.channelLabel{2:end} '.png']));
end

%% radial probability of being positive channel by channel (MP Only)

markerChannels = 2:4;
for channel = markerChannels
    figure()
    colors = lines(numel(manMeta.conditions));

    lw = 2;
    fs = 20;
    legendstr = [];
    p = [];
    for condi = 1:numel(meta.conditions)
        Pspe = Pall{condi}(channel,:);
        Pstd_spe = Pstdall{condi}(channel,:);
        xspe = xall{condi};

        ptemp = plot(xspe,Pspe,'LineWidth',lw,'Color',colors(condi,:));
        legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
        p = [p,ptemp];
        hold on
        good = ~isnan(Pspe);
        fill([xspe(good),fliplr(xspe(good))],...
            [Pspe(good) + Pstd_spe(good), fliplr(Pspe(good) - Pstd_spe(good))],...
            colors(condi,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    xlim([0 max(xall{condi})]);
    %ylim([0 1]);
    xlabel('edge distance ( um )')
    ylabel('positive fraction ');
    legend(p,legendstr,'FontSize',15,'Location','best')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    saveas(gcf, fullfile(dataDir, ['radialpositive_combine_', meta.channelLabel{channel} '.png']));
end

%% count populations

combo = [4 3 2];
conditionsidx = 1:numel(manMeta.conditions);
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% make pretty scatter plot

for condi = 1:numel(meta.conditions)

    close all;
    options = struct();

    options.selectDisType = [1, 1, 1, 1]; % 1 - nuc, 2 - cyto, 3 - NCratio, 4 - CNratio
    options.channelThresholds = stats.thresholds; % re-define channelThres when plot cyto or NCratio
    %options.checkPos = badidxCon; % whether you want to check defined cells in scatter plot
    options.showThreshold = [false false false false];
    options.minimal = true;
    options.edgeDistance = false;

    options.imageType = imageType;
    options.conditionIdx = condi;
    options.channelCombos = {[4 2 3],[4 3 2],[3 2 4]};%, [4 2 3]};%, [4 2 3]}; % {[3 2 4], [4 3 2], [4 2 3]}

    options.axisLabel = meta.channelLabel;
    options.channelMax = exp([3 3 3 3])-1; %xlim or ylim
    options.log1p = [true true true true];
    %options.log1p = [false false false false];
    if strcmp(imageType,'MP')
        options.radiusMicron = max(xall{condi}); %positions(meta.conditionStartPos(condi)).radiusMicron;
    end
    options.conditionsCombined = false;
    %options.showThreshold = [false false true false];

    results = scatterMicropattern(stats, meta, dataDir, options);
    %writetable(results.dataTable, fullfile(dataDir, ['scatterData_' meta.conditions{condi} '.csv']))
end


%% pretty channel combo pie (MP Only)

options = struct();
options.channels = [4 3 2];
options.pieOrder = [1 2 3];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.01 0.995; 0.01 0.995; 0.01 0.995];
options.tol = [0.1 0.995; 0.1 0.995; 0.1 0.995; 0.1 0.995];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.scalebar = false;
options.color = 'RGB';

posidx = 14:-1:1; % the first position sets the lookup table

for pi = posidx

    condi = find(pi >= meta.conditionStartPos,1,'last');
    cpi = pi - meta.conditionStartPos(condi) + 1;
    %fname = filelist{condi}(cpi).name;
    %subDir = filelist{condi}(cpi).folder;
    %subDir = fullfile(dataDir, 'MIP');
    
    if pi == posidx(1)
        Ilim = micropatternPieVis(dataDir, positions(pi), options);
    else
        options.Ilim = Ilim;
        micropatternPieVis(dataDir, positions(pi), options);
    end
end

%% ANALYZE SCRATCHES

for pi = 1:14

    % load mask
    rawmask = load(sprintf('stitched_p%.4d_w0000_t0000_masks.mat',pi-1));
    rawmask = ~min(rawmask.bgmask,[],3);

    % clean up
    mask = imclose(rawmask,strel('disk',3));
    mask = bwareaopen(mask,3000); % remove small objects
    mask = imclose(mask,strel('disk',10));
    mask = ~bwareaopen(~mask,10000); % fill small holes
    %imshow(cat(3, uint8(mask), rawmask, mask)*255)

    % get largest connected components of difference between mask and convex
    % hull

    %imshow(cat(3, uint8(mask), rawmask, bwconvhull(mask))*255)

    scratchmask = imopen(bwconvhull(mask) - mask,strel('disk',3));
    CC = bwconncomp(scratchmask);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [sortnumPixels,I] = sort(numPixels,'descend');
    scratchmask = false(size(scratchmask));
    scratchmask(cat(1,CC.PixelIdxList{I(sortnumPixels > 40000)})) = true;
    %imshow(scratchmask);

    % exclude anything close to the edge
    radialMargin = -30;
    position = positions(pi);
    xres = position.radiusMicron / position.radiusPixel;
    center = position.center;
    Rmax = uint16(position.radiusPixel + radialMargin/xres); % micron margin
    [X,Y] = meshgrid(1:size(mask,2),1:size(mask,1));
    R = sqrt((X - center(1)).^2 + (Y - center(2)).^2);
    F = atan2((X - center(1)),-(Y - center(2)));
    disk = R < Rmax;
    %imshow(cat(3,uint8(mask),disk,mask)*255)

    scratchmask = scratchmask.*disk;
    %imshow(scratchmask)

    % distance from scratch
    D = bwdist(scratchmask).*mask;
    %imshow(D,[])

    % load ERK image
    %------------------------

    % self stitched files start with stitched
    if strcmp(position.filenameFormat(1:8), 'stitched') ||...
        strcmp(position.filenameFormat(1:10), 'FFstitched') 
        % one version of the code (not current) used FF as a prefix for
        % flatfield corrected

        s = strsplit(position.filename,{'_'});

        if ~isempty(dir(fullfile(dataDir,'filtered',[s{1} '_filtered_MIP_' s{2} '*'])))
            disp('using filtered MIP')
            prefix = [s{1} '_filtered_MIP_' s{2}];
            subDir = 'filtered';
        else
            prefix = [s{1} '_MIP_' s{2}];
            subDir = 'MIP';
        end
        fnameFormat = [prefix '_w%.4d_t0000.tif'];

    % fusionstitcher case: [fname]_FusionStitcher.ims
    else
        s = strsplit(position.filename,{'_FusionStitcher','.ims'});

        if ~isempty(dir(fullfile(dataDir,'filtered',[s{1} '_filtered_MIP*'])))
            disp('using filtered MIP')
            prefix = [s{1} '_filtered_MIP'];
            subDir = 'filtered';
        else
            prefix = [s{1} '_MIP'];
            subDir = [prefix '_zslices'];
        end
        fnameFormat = [prefix '_w%.4d.tif'];
    end

    imgs = {};
    Ilim = {};
    for cii = 4

        fullfname = fullfile(dataDir, subDir, sprintf(fnameFormat, cii-1));
        try
            img = imread(fullfname);
        catch
            img = imread([fullfname(1:end-4) '.jpg']);
        end
    end

    %imshow(imadjust(img),[])

    smax = 1200; %round(max(D(:))

    nbins = 50;
    ds = smax/nbins;
    X = zeros([nbins 2]);

    for i = 0:nbins-1
        X(i+1,2) = mean(img((D > ds*i) & (D < ds*(i+1))));
        X(i+1,1) = mean(D((D > ds*i) & (D < ds*(i+1))));
    end

    if sum(scratchmask(:)) > 10^4
        graphs{pi} = X;
    else
        graphs{pi} = [];
    end
end

%figure, plot(x,y)

%% visualize

figure, 
meany = {};
stdy = {};
ymax = 0;
ymin = 10^6;

conditionindices = [3 1];
colors = lines(numel(conditionindices));

for condi = conditionindices
    
    ally = {};

    for pi = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1

        if ~isempty(graphs{pi})
            %plot(graphs{pi}(:,1), graphs{pi}(:,2))
            ally{pi} =  graphs{pi}(:,2);
        end
    end
    
    ally = cat(2,ally{:});
    disp(size(ally,2));
    meany{condi} = nanmean(ally,2);
    
    ymin = min(ymin, min(meany{condi}));
    ymax = max(ymax, max(meany{condi}));
    stdy{condi} = nanstd(ally,0,2);
end

%figure,
lw = 3;
fs = 26;
clf
hold on
x = graphs{1}(:,1)'*meta.xres;
i = 0;
for condi = conditionindices
    i = i+1;
    normmean = (meany{condi}-ymin)/(ymax-ymin);
    plot(x, normmean, 'LineWidth',lw,'Color',colors(i,:))
end
i = 0;
for condi = conditionindices
    i = i+1;

    normmean = (meany{condi}-ymin)/(ymax-ymin);
    normstd = stdy{condi}/(ymax-ymin);
    
    plerr = normmean + normstd;
    minerr = normmean - normstd;
    good = ~isnan(plerr);
    fill([x(good),fliplr(x(good))],[plerr(good)', fliplr(minerr(good)')],...
        colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
end
hold off
cleanSubplot(fs,lw);
xlim([0 300]);
ylim([-0.1 1.1]);
yticks([0 0.5 1]);
xlabel('distance from scratch (\mum)');
ylabel('pERK intensity (a.u.)   ');
axis square;
legend({'scratch','+FGFRi'})

saveas(gcf, ['pERKprofile_scratch' '.png'])
%% Export scratch assay pERK profile data to CSV files

% Create output directory for CSV files
csvDir = fullfile(dataDir, 'scratch_profile_csv');
if ~exist(csvDir, 'dir')
    mkdir(csvDir);
end

conditionindices = [3 1];
conditionNames = {'scratch', '+FGFRi'};

% Collect all data for combined export
allData = [];

for i = 1:numel(conditionindices)
    condi = conditionindices(i);
    
    % Get x values (distance)
    x = graphs{1}(:,1)' * meta.xres;
    x = x(:);  % Force column vector
    
    % Get mean and std
    normmean = (meany{condi} - ymin) / (ymax - ymin);
    normstd = stdy{condi} / (ymax - ymin);
    
    % Also save raw (unnormalized) values
    raw_mean = meany{condi};
    raw_std = stdy{condi};
    
    % Remove NaN values
    good = ~isnan(normmean);
    
    % Create table for this condition
    conditionName = conditionNames{i};
    
    tempTable = table(...
        x(good), ...
        normmean(good), ...
        normstd(good), ...
        raw_mean(good), ...
        raw_std(good), ...
        repmat({conditionName}, sum(good), 1), ...
        'VariableNames', {'Distance_um', 'Normalized_Intensity', 'Normalized_Std', 'Raw_Intensity', 'Raw_Std', 'Condition'});
    
    % Append to combined data
    allData = [allData; tempTable];
end

% Save combined CSV
csvFilename = fullfile(csvDir, 'pERK_profile_scratch_all_conditions.csv');
writetable(allData, csvFilename);
fprintf('Saved combined file: %s\n', csvFilename);

%%
pi = 7;
options = struct();
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.color = 'RGB';
options.channels = [4 3];
options.pieOrder = [1 2];
micropatternPieVis(dataDir, positions(pi), options);

%% multi condition pie (MP Only)

options = struct();
options.channels = [4 3];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.1 0.99; 0.1 0.99; 0.6 0.95];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.color = 'RGB';
options.scalebar = true;
options.overlayweights = [0.9 0.7];

options.positionIdx = [8 2];

micropatternPieVisConditions(dataDir, positions, options, meta);


%%
% ------------------------------------------------------------------------
% cell density and confluence
% ------------------------------------------------------------------------

%% density plot condition by condition (MP Only)

rAll = cell(numel(meta.conditions),1);
densityAll = cell(numel(meta.conditions),1);

numberAll = zeros(manMeta.posPerCondition(1), numel(manMeta.conditions));
for condi = 1:numel(meta.conditions)
    rAll{condi} = res{condi}.r';
    density = [];
    for pos = 1:manMeta.posPerCondition
        posidx = (condi-1)*manMeta.posPerCondition(1)+pos;
        densityTem = res{condi}.celldensity{posidx}';
        density = [density,densityTem];
        numberAll(pos,condi) = positions(posidx).ncells;
    end
    densityAll{condi} = density;
end

PI = 3.1415926;
densityColonyAll = numberAll/(PI*350^2)*10^6;
numberAllMean = mean(numberAll,1);
numberAllStd = std(numberAll,0,1);

c = lines(5);
fs = 20;
p = [];
figure()
for condi = 1:numel(meta.conditions)
    r = rAll{condi};
    densityMEAN = mean(densityAll{condi}*10^6,2); % cells/mm^2
    densitySTD = std(densityAll{condi}*10^6,0,2); % cells/mm^2
    hold on
    errorbar(r, densityMEAN, densitySTD,'--','LineWidth',0.5, 'Color', c(condi,:));
    ptemp = plot(r, densityMEAN,'-','LineWidth',3, 'Color', c(condi,:), 'DisplayName', manMeta.conditions{condi});
    hold off
    p = [p,ptemp];
end
xlim([0 350])
ylim([0 30000])
legend(p,'Location','best');
xlabel('edge distance ( um )')
ylabel('cell density (cells/mm^2)');
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);
saveas(gcf, ['radialProfile_density' '.png'])

%% Calculate density & confluency of not confluent position

posNoConfluent = [1,5,9];
Area = zeros(numel(positions),1);
Confluency = ones(numel(positions),1);
Density = zeros(numel(positions),1);
for pi = posNoConfluent
    condi = ceil(pi / meta.posPerCondition);
    cpi = pi - (condi-1)*meta.posPerCondition;
    img = max(positions(pi).loadImage(dataDir, 0, 1),[],3);
    img = imadjust(img,stitchedlim(img));
    imgBW = imbinarize(img);
    imgBWCLEAN = bwareaopen(imgBW,1000); % remove small particles
    diskSize = 30; % play with disk size
    imgBW2 = imfill(imclose(imgBWCLEAN,strel('disk',diskSize)),'holes');
    figure()
    imshowpair(imgBWCLEAN,imgBW2,'montage')
    saveas(gcf,['Area_' manMeta.conditions{condi} '_Position_' num2str(pi) '.jpg']);
    Area(pi) = bwarea(imgBW2)*meta.xres*meta.yres/10^6;
end
for pi = 1:numel(positions)
    if Area(pi) == 0
        Area(pi) = stats.area*100;
    else
        Confluency(pi) = Area(pi)/(stats.area*100);
    end
    Density(pi) = positions(pi).ncells/Area(pi); % cells/mm^2
end
DensityAll = zeros(numel(manMeta.conditions),1);
DensitySTD = zeros(numel(manMeta.conditions),1);
coloniesPerSize = cellfun(@numel,filelist);
posbins = [0 cumsum(coloniesPerSize)];
for condi = 1:numel(manMeta.conditions)
    DensityAll(condi) = mean(Density(posbins(condi)+1:posbins(condi+1)));
    DensitySTD(condi) = std(Density(posbins(condi)+1:posbins(condi+1)));
end
