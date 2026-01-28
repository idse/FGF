% VERSION HISTORY
%
% 231115 : Idse - change to deal with self-stitched files
%

clear; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);

% define channels and conditions
%--------------------------------------
% define micropatterned 'MP' or 'Disordered'
imageType = 'MP'; 
% margin in microns to keep around the nominal radius (for masking cells)
Rmarg=55; 

% different cases for input files
%--------------------------------------

% self-stitched - channels are split
if ~isempty(dir(fullfile(dataDir,'stitched*.tif')))
    filelistAll = {};
    filenameFormat = "stitched_p%.4d_w%.4d_t%.4d.tif";
    xsectfname_prefix_format = 'stitched_p%.4d';
    
% Dragonfly stitched - channels are not split
elseif ~isempty(dir(fullfile(dataDir,'*FusionStitcher*.ims')))
    filelistAll = dir(fullfile(dataDir,'*FusionStitcher*.ims'));
    % FINISH
end

% define files per condition
%--------------------------------------

filelist = {};

% case 1: fixed number of positions per condtion 
%posPerCondition = 4;
% manMeta.posPerCondition = repmat(posPerCondition, [1 nWells]);
% manMeta.nPositions = posPerCondition*nWells;
% conditionStartPos = 1 + [0 cumsum(manMeta.posPerCondition)];
% conditionStartPos = conditionStartPos(1:end-1);

if ~isempty(filelistAll)
    for i = 1:numel(manMeta.conditions)
        filelist{i} = filelistAll((i-1)*posPerCondition+1:i*posPerCondition);
    end
% else
%     filelist = {};
% 
%     for condi = manMeta.nWells
%        filelist{condi} = []; 
%     end
% 
%     for pi = 1:manMeta.nPositions
%         condi = find(pi >= conditionStartPos,1,'last');
%         filelist{condi} = [filelist{condi} sprintf(filenameFormat, pi)];
%     end
end

% % case 2: define files per condition manually 
% % required if variable number of positions per conditions
% filelist{1} = dir(fullfile(dataDir,'*80K*.ims'));
% filelist{2} = dir(fullfile(dataDir,'*160K*.ims'));
% filelist{3} = dir(fullfile(dataDir,'*240K*.ims'));
% filelistAll = cat(2,filelist(:));
% 
% % define number of positions for each condition
% manMeta.posPerCondition = cellfun(@numel,filelist);

% starting position of each conditioni is meta.conditionStartPos

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

% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = ones([1 numel(meta.conditions)])*700/2; 


% ratio between z and xy resolution
scalefactor = meta.zres / meta.xres;


%% cleanup stack for specific positions

zcorrection = true;
zcorrfactor = [1.02 1.02 1 1];
mediansize = [5 5];

tic
%parpool(3)
%posidx = 5;%[1,5,9];
posidx = [1];
for i = 1:numel(posidx)
%for i = 1:numel(posidx)
    fname = char(sprintf(filenameFormat, posidx(i)-1));
    cleanupStack(dataDir, fname, zcorrection, zcorrfactor, mediansize)
end
toc


%micropatternPieVisConditions(dataDir, positions, options, meta);

%% make xsections of different channels

pList = [3,5,9];
RGBset = [2 3 4];

% parameters for cross-section
tol = [0.01 0.995; 0.01 0.995; 0.01 0.995; 0.01 0.995];
yinds = {1350:1360, 1300:1310, 1275:1285, 1250:1260, 1200:1210, 1150:1160, 1100:1110, 1050:1060,1025:1035, 1000:1010, 975:985, 950:960, 925:935, 900:910, 875:885}; %change yinds if preferred


pList = 9;
yinds = {1350:1360};
RGBset = [3];

count = 0; % use first position to define lim
for pi = pList
    count = count + 1;
    % self-stitched (channels split)
    if filenameFormat.startsWith('stitched')
        i = 0;
        for ci = 0:meta.nChannels-1
            i = i + 1;
            fullfname = fullfile(dataDir,sprintf(filenameFormat, pi-1, ci, 0));
            [fpath,barefname,ext] = fileparts(fullfname);
            barefname = char(barefname);
            filteredfname = fullfile(fpath, 'filtered', ['stitched_filtered' barefname(9:end) char(ext)]);
            
            if exist(filteredfname)
                disp('using filtered stack');
                fullfname = filteredfname;
            else
                warning('no filtered stack');
            end
            
            stack{i} = readStack(fullfname);
            % readStack mistakes z for t in stitched tiffs
            stack{i} = permute(stack{i}, [1 2 3 5 4]); 
        end
        img = cat(3,stack{:});
        xsectfname_prefix = sprintf(xsectfname_prefix_format, pi-1);
        
    % not self-stitched 
    else
        fname = filelistAll{pi}.name;
        img = readStack(fullfile(dataDir,fname));
        xsectfname_prefix = fname(1:end-4);
    end

    for di = 1:2
        if di == 1
            direction = 'x';
        else
            direction = 'y';
        end
        for i = 1:numel(yinds)

            options = struct('tol',tol, 'index', yinds{i}, 'RGBset', RGBset,...
                'color', 'GM','direction',direction);

            if count == 1
                disp('setting Ilim');
                options.Ilim = makeXsection(img, xsectfname_prefix, meta, options);
            else
                makeXsection(img, xsectfname_prefix, meta, options);
            end
            close;
        end
    end
end


%--------------------------------------------------------------------------
% INITIAL QC
%--------------------------------------------------------------------------

%% overlay quantification for QC on XY MIP

piSerious = 5;%[1,5,9];
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
    c = z;
    
    z = positions(pi).cellData.Z;
    z = (z-min(z))/(max(z)-min(z));
    
    figure()
    imshow(imgGREY)
    hold on
    %s = scatter(XY(:,1), XY(:,2), 50, c*[1, 1, 1], "filled");
    s = scatter(XY(:,1), XY(:,2), 50, z*[1, 1, 1], "filled");
    hold off
    %s.MarkerEdgeColor = 'm';
    title(['Intensity in ',meta.channelLabel{ci}])
    saveas(gcf,['Intensity_check_' meta.channelLabel{ci} '_' meta.conditions{condi} '_Pos' num2str(pi) '.jpg']);
end

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


%% get stats

stats = cellStats(positions, meta, positions(1).dataChannels);
% automaticly get thresholds
confidence = 0.95;
conditions = 1:numel(meta.conditions);
whichthreshold = []; %[1 1 2]
stats.getThresholds(confidence, conditions, whichthreshold);

%%
% Make nucHistogram channel by channel
conditionIdx = 1:numel(meta.conditions); % [1,3,5] Specify conditions for Histograms
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
for channelIdx = 1:numel(meta.channelLabel)
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
%%
save(statsfile, 'stats');
%stats.exportCSV(meta); export to CSV

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

%% make intensity radial profile (MP Only)

% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', method 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD


% options = struct('nucChannels', [1 2 3], ...'ratioMode', 'N:C',...
%     'normalize', true, 'std', true, 'FontSize', 15, 'legend',true,...
%     'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
%     'log1p',false);

options = struct('nucChannels', [2], ...
    'cytChannels', [1],...'ratioMode', 'N:C',...
    'normalize', true, 'std', true, 'FontSize', 15, 'legend',false,...
    'pointsPerBin', 300, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
    'log1p',false);

% options = struct('ratioChannels', 1, 'ratioMode', 'N:C',...
%     'normalize', true, 'std', true, 'FontSize', 15, 'legend',true,...
%     'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
%     'log1p',false);

res = {};
% first round without normalization (to collect overall max and min of all), keep normalize = false



conditionindices = [1 2 3];

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
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    figure()
    res{condi} = plotRadialProfiles(stats, meta, options);
    %ylim([0 1.1]);
    %ylim([-0.0195 -0.016]);
    xlim([0 radii(condi)]);
    axis square
    title(meta.conditions{condi})
    saveas(gcf, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.png'])
    %writetable(res{condi}.datatable, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.csv'])
    close;
end

% radial intensity - conditions combined

colors = lines(numel(conditionindices));
lw = 3;
fs = 32;
    
markerChannels = options.nucChannels+1;
for channel = markerChannels
    figure()
    legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.nuc_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        legendstr = [legendstr,meta.conditions(condi)];
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
    legendstr = {'BMP4 only','+FGFRi','+FGFRi+dox'};
    %legend(p,legendstr,'FontSize',20,'Location','best')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    %title(meta.channelLabel{channel})
    axis square;
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', legendstr{:}, '_' meta.channelLabel{channel} '.png']));
end

% CYTOPLASMIC CHANNELS
markerChannels = options.cytChannels+1;
for channel = markerChannels
    figure()
    legendstr = [];
    p = [];
    
    clf
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.cyt_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        legendstr = [legendstr,meta.conditions(condi)];
    end
    i = 0;
    for condi = conditionindices
        i = i+1;
        
        plerr = res{condi}.cyt_profile(:,channel) + res{condi}.cyt_profile_colstd(:,channel);
        minerr = res{condi}.cyt_profile(:,channel) - res{condi}.cyt_profile_colstd(:,channel);
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
    legend(p,legendstr,'FontSize',15,'Location','best')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    title(meta.channelLabel{channel})
    axis square;
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', legendstr{:}, '_' meta.channelLabel{channel} '.png']));
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
posList = [1,5,9];
norMethodPos = 1; % 1 - normalize individually; 2 - normalize according to all position
options = struct('nucChannels', 1:3, 'ratioMode', 'C:N',...
    'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 100, 'interpmethod', 'nearest', 'stdMethod', 'perColony');
% first round without normalization (to collect overall max and min of all), keep normalize = false
count = 0;
for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition(1));
    resOP = {};
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    metaOP.nWells = 1;
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
    metaOP.nWells = 1;
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

condi = 2; % condition index
cpi = 1; % position within condition index
s = 30;

pi = cpi + meta.conditionStartPos(condi);
%subDir = filelist{condi}(pi).folder;

%channelThresholds = [0, stats.thresholds(2), stats.thresholds(3), stats.thresholds(4)]; % use automatic thresholds
channelThresholds = [0 500 250 400]; % manually adjust thresholds

nucLevel = positions(pi).cellData.nucLevel;
background = positions(pi).cellData.background;
positive = {};
for i = 1:meta.nChannels
    positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
end

for ci = 3%:numel(manMeta.channelLabel)
    
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


%%
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


markerChannels = 3;
for channel = markerChannels
    figure()
    colors = lines(numel(meta.conditions));

    lw = 3;
    fs = 32;
    legendstr = [];
    p = [];
    for condi = [2 1 3]
        Pspe = Pall{condi}(channel,:);
        Pstd_spe = Pstdall{condi}(channel,:);
        xspe = mean([res{condi}.trueradius{:}]) - xall{condi};

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
    ylim([0 1.1]);
    xlabel('edge distance (\mum)')
    ylabel('positive fraction ');
    legendstr = {'BMP4','+FGFRi','+dox'};
    %   title(meta.channelLabel{channel});
    legend(p,legendstr,'FontSize',20,'Location','NorthEast')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', lw);
    axis square
    box off
    saveas(gcf, fullfile(dataDir, ['radialpositive_combine_', meta.channelLabel{channel} '.png']));
end
%% Export radial probability data to CSV files

% Create output directory for CSV files
csvDir = fullfile(dataDir, 'radial_data_csv');
if ~exist(csvDir, 'dir')
    mkdir(csvDir);
end

markerChannels = 3;
for channel = markerChannels
    
    % Initialize combined data structure
    allData = [];
    
    for condi = [2 1 3]
        % Extract data
        Pspe = Pall{condi}(channel,:);
        Pstd_spe = Pstdall{condi}(channel,:);
        xspe = mean([res{condi}.trueradius{:}]) - xall{condi};
        
        % Remove NaN values
        good = ~isnan(Pspe);
        
        % Create table for this condition
        conditionName = meta.conditions{condi};
        channelName = meta.channelLabel{channel};
        
        tempTable = table(...
            xspe(good)', ...
            Pspe(good)', ...
            Pstd_spe(good)', ...
            repmat({conditionName}, sum(good), 1), ...
            repmat({channelName}, sum(good), 1), ...
            'VariableNames', {'Distance_um', 'Positive_Fraction', 'Std_Dev', 'Condition', 'Channel'});
        
        % Append to combined data
        allData = [allData; tempTable];
    end
    
    % Save combined CSV for this channel
    csvFilename = fullfile(csvDir, ['radial_probability_', meta.channelLabel{channel}, '_all_conditions.csv']);
    writetable(allData, csvFilename);
    fprintf('Saved: %s\n', csvFilename);
    
end

%% count populations

combo = [4 3 2];
conditionsidx = 1:numel(meta.conditions);
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% cell number statistical significance

cell_counts = counts.Ncells([5:8, 1:4, 9:12]); %reshape(counts.Ncells, [4 3]);

groups = [repmat({'control'}, 1, 4), ...
            repmat({'Fi'}, 1, 4), ...
          repmat({'Fid'}, 1, 4)];

% Step 1: One-way ANOVA
[p, tbl, stats] = anova1(cell_counts, groups, 'off');  % 'off' suppresses plot
disp(['ANOVA p-value = ', num2str(p)])

% Step 2: Tukey's HSD post-hoc test
results = multcompare(stats);
close;

% Find which comparisons are statistically significant
significant = results(:,6) < 0.05;
disp('Significant comparisons:')
disp(results(significant, :))

% Extract comparison info
% multcompare output columns:
% [group1, group2, lowerCI, meanDiff, upperCI, p-value]

sig_pairs = {};
sig_pvals = [];

for i = 1:size(results, 1)
    g1 = results(i, 1);
    g2 = results(i, 2);
    pval = results(i, 6);

    % Store as sigstar expects: {[x1 x2], ...}, [p1, p2, ...]
    sig_pairs{end+1} = [g1 g2];
    sig_pvals(end+1) = pval;
end

%% cell number bar plot

counts.NcellsTotal([2 1 3]);
fs = 38;
lw = 3;
fnameprefix = 'cellnumbers_bar';
vals = counts.NcellsTotal([2 1 3])'./meta.posPerCondition([2 1 3]);
errs = counts.NcellsStd(conditionsidx)';
ylabelstr = 'cells / colony';

labels = {'BMP4','+FGFRi','+FGFRi+dox'};
%labels = {'x10^3','+Fi','+dox'};

errorbar_groups(vals, errs,...
            'bar_names', labels,...
            'bar_width',0.75,'errorbar_width',0.5,...
            'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',lw});
ylim([0 max(vals(:))+max(errs(:))]);

sig_pvals(sig_pvals > 0.05) = NaN;
sigstar(sig_pairs, sig_pvals);
% make the stars bigger
text = findall(gca, 'Type', 'Text');
for s = text'
    txt = s.String;
    if ischar(txt) && (contains(txt, '*') || contains(txt, 'p =') || contains(txt, 'n.s.'))
        s.FontSize = 30;  
    end
end


bgc='w';
fgc='k';
ylabel(ylabelstr);
set(gcf,'color',bgc);
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',bgc);
yt = [0 1 2 3 4 5]*1000; % get(gca,'YTick');
set(gca,'YTick',yt)
set(gca,'YTickLabel', yt/1000)
axis square


saveas(gcf, fullfile(dataDir, [fnameprefix num2str(conditionsidx) '.png'])); 
%% Export all bar plot data to a single CSV

% Create a cell array to store all data with clear section headers
export_cell = {};

% ===== SECTION 1: ANOVA Results =====
export_cell{end+1, 1} = 'ANOVA RESULTS';
export_cell{end+1, 1} = 'Test';
export_cell{end, 2} = 'P-Value';
export_cell{end, 3} = 'Significant (p<0.05)';

export_cell{end+1, 1} = 'One-way ANOVA';
export_cell{end, 2} = p;
export_cell{end, 3} = p < 0.05;

% Add blank row
export_cell{end+1, 1} = '';

% ===== SECTION 2: Summary Statistics =====
export_cell{end+1, 1} = 'SUMMARY STATISTICS';
export_cell{end+1, 1} = 'Condition';
export_cell{end, 2} = 'Mean (Cells/Colony)';
export_cell{end, 3} = 'Std Error';
if exist('meta', 'var') && isfield(meta, 'posPerCondition')
    export_cell{end, 4} = 'N (Colonies)';
end

for i = 1:length(labels)
    export_cell{end+1, 1} = labels{i};
    export_cell{end, 2} = vals(i);
    export_cell{end, 3} = errs(i);
    if exist('meta', 'var') && isfield(meta, 'posPerCondition')
        export_cell{end, 4} = meta.posPerCondition([2 1 3]); % Match the order
    end
end

% Add blank row
export_cell{end+1, 1} = '';

% ===== SECTION 3: Pairwise Comparisons =====
export_cell{end+1, 1} = 'PAIRWISE COMPARISONS (Tukey HSD)';
export_cell{end+1, 1} = 'Comparison';
export_cell{end, 2} = 'Mean Difference';
export_cell{end, 3} = '95% CI Lower';
export_cell{end, 4} = '95% CI Upper';
export_cell{end, 5} = 'P-Value';
export_cell{end, 6} = 'Significance';

% Map group numbers to condition names
group_map = containers.Map([1 2 3], labels);

for i = 1:size(results, 1)
    g1 = results(i, 1);
    g2 = results(i, 2);
    pval = results(i, 6);
    
    export_cell{end+1, 1} = sprintf('%s vs %s', group_map(g1), group_map(g2));
    export_cell{end, 2} = results(i, 4);  % Mean difference
    export_cell{end, 3} = results(i, 3);  % Lower CI
    export_cell{end, 4} = results(i, 5);  % Upper CI
    export_cell{end, 5} = pval;
    
    % Add significance stars
    if pval <= 0.0001
        export_cell{end, 6} = '****';
    elseif pval <= 0.001
        export_cell{end, 6} = '***';
    elseif pval <= 0.01
        export_cell{end, 6} = '**';
    elseif pval <= 0.05
        export_cell{end, 6} = '*';
    else
        export_cell{end, 6} = 'n.s.';
    end
end

% Convert cell array to table and export
export_table = cell2table(export_cell);

% Write to CSV
csv_filename = fullfile(dataDir, [fnameprefix num2str(conditionsidx) '_complete_analysis.csv']);
writetable(export_table, csv_filename, 'WriteVariableNames', false);

disp(['Exported complete analysis to: ' csv_filename]);

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
options.channels = [2 3 4];
options.pieOrder = [1 2 3];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.01 0.995; 0.01 0.995; 0.01 0.995];
options.tol = [0.1 0.995; 0.1 0.995; 0.1 0.995; 0.1 0.995];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.scalebar = false;
options.color = 'CMY';

posidx = [10]%1:12; % the first position sets the lookup table

% 3 5 10

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


%%
pi = 1;
options = struct();
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.color = 'MG';
options.channels = [2 4];
options.pieOrder = [1 2];
micropatternPieVis(dataDir, positions(pi), options);


%% multi condition pie (MP Only)

options = struct();
options.channels = [2 3 4];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.1 0.8; 0.1 0.9; 0.01 0.99];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.color = 'RGB';
options.scalebar = true;
options.overlayweights = [0.9 0.7]; % weight of colors vs DAPI image 

% options.positionIdx = [1 5 9];
% options.positionIdx = [3 6 11];

for i = 0%:3
    %options.positionIdx = meta.conditionStartPos([2,4,5])+i; % 14h
    %options.positionIdx = meta.conditionStartPos([2,7,8])+i; % 24h
    %options.positionIdx = meta.conditionStartPos([2,9,10])+i; % 2h
    %options.positionIdx = meta.conditionStartPos([2 3 9 10])+i; % 2h + 10min
    %options.positionIdx = [26 29 8 9];
    options.positionIdx = [5 10];
    
    micropatternPieVisConditions(dataDir, positions, options, meta);
end


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
