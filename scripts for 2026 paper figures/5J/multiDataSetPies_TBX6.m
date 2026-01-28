clear; close all; clc

%%

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
[parentDir,~,~] = fileparts(scriptPath);
[grandparentDir,~,~] = fileparts(parentDir);
[greatgrandparentDir,~,~] = fileparts(grandparentDir);

baseDirs = fullfile(scriptPath,...
    {'20250620_PGP1 shRNA Scr F17 1 6_SOX2 TBX6 TBXT',...
    '20240812_ESI shRNA Scr FGF4_SOX17 TBX6 BRA'});

expnums = [7 10];

nb = length(baseDirs);

subdirs = cell(1,nb);
channelLabels = cell(1,nb);
channelDirs = cell(1,nb);
channelIds = cell(1,nb);
pos = cell(1,nb);
for bi = 1:nb
    baseDir = baseDirs{bi};
    listing = dir(fullfile(baseDir));
    listing = listing([listing.isdir]);
    subdirs{bi} = {listing.name};

    load(fullfile(baseDir,'meta.mat'),'meta')
    load(fullfile(baseDir,'positions.mat'),'positions')
    channelLabels{bi} = renameDuplicateChannels(meta.channelLabel);
    channelDirs{bi} = cumsum(contains(meta.channelLabel,'DAPI'));
    pos{bi} = positions;

    dirids = unique(channelDirs{bi});
    IDs = NaN(size(channelDirs{bi}));
    for di = 1:length(dirids)
        mask = channelDirs{bi} == dirids(di);
        cs = cumsum(mask);
        IDs(mask) = cs(mask);
    end
    channelIds{bi} = IDs;
end

chans = intersect(channelLabels{:});
chans = chans(~(contains(chans,'DAPI') & ~strcmp(chans,'DAPI_1')));
nchan = length(chans);


%% look at individual images
bare = 'stitched_MIP_p%.4d_w%.4d_t0000.jpg';
chan = 'TBXT';
tol = [0.2 0.995];
scale = 1;
lims = NaN(nb,2);

pis = {[1 4],[4 5]};
nrow = length(pis); ncolumn = max(cellfun(@length,pis));
titles = {{'ctrl shRNA','FGF17 shRNA'},{'ctrl shRNA','FGF4 shRNA'}};
figure('Position',figurePosition(1350,1200))
for bi = 1:nb
    colis = pis{bi};
    %cii = find(strcmp(channelLabels{bi},chan));
    cii = 4;
    di = channelDirs{bi}(cii);
    chanid = channelIds{bi}(cii);
    for ii = 1:length(colis)
        coli = colis(ii);
	    fname = fullfile(baseDirs{bi},subdirs{bi}{di},'MIP',sprintf(bare,coli-1,chanid-1));
        im = imread(fname);
        if ii == 1
            lim = stitchedlim(im,tol);
            lim(2) = lim(1) + scale*(lim(2) - lim(1));
            % lim(2) = scale*lim(2);
            lims(bi,:) = lim;
        end
        subplot_tight(nrow,ncolumn,(bi - 1)*ncolumn + ii,0.03)
        imshow(imadjust(im,lims(bi,:)))
        title(titles{bi}{ii},'FontWeight','normal')
    end
end


%%
savedir = fullfile(baseDirs{1},'pieSlices');
if ~exist(savedir,'dir'), mkdir(savedir); end

refcols = [2 2];
%piecols = [2 8 10];
piecols = [3 4 5];
piedirs = [1 1 2];
N = length(piecols);

%set individual tolerances
tols = repmat([0.01 0.99],nchan,1);

tols(strcmp(chans,'TBXT'),:) = [0.2, 0.995];


%scale the limit set based on tolerances for pAKT to get high values based
%on the control without trying to use an unreasonable percentage for the
%tolerance in the control condition
scales = ones(nchan,1);

saveprev = sprintf('pieSlices');

PI = atan(1)*4;
margin = 200;
radialMargin = 20;
outsidemargin = 150;
padsize = round(outsidemargin/2);

bare = 'stitched_filtered_MIP_p%.4d_w%.4d_t0000.tif';

nameidxs = [expnums(piedirs); piecols-1];
pref = ['pieSlices',sprintf('_ExpX%.1d_p%.4d',nameidxs(:))];

%%
bi = 1;
coli = 1;

xres = pos{bi}(coli).radiusMicron/pos{bi}(coli).radiusPixel;
Rmax = round(pos{bi}(coli).radiusPixel + radialMargin/xres); % micron margin
Rcrop = Rmax + padsize;
cropsize = 2*Rcrop + 1;

[X,Y] = meshgrid(1:cropsize,1:cropsize);
R = sqrt((X - (Rcrop+1)).^2 + (Y - (Rcrop+1)).^2);
F = atan2((X - (Rcrop+1)),-(Y - (Rcrop+1)));
disk = R > Rmax;

mask = cell(1,N);
for ii = 1:N
    mask{ii} = F < -pi + (ii-1)*2*pi/N | F > -pi + ii*2*pi/N;
    mask{ii} = imdilate(mask{ii},strel('disk',10));
end
lines = mask{1};
for ii = 2:N
    lines = lines & mask{ii};
end
%%
close all
figure
for ci = 1
    scale = scales(ci);
    chan = 'TBXT';
    lims = NaN(nb,2);
    for bi = 1:nb
        coli = refcols(bi);
        %cii = find(strcmp(channelLabels{bi},chan));
        cii = 4;
        di = channelDirs{bi}(cii);
        chanid = channelIds{bi}(cii);
        fname = fullfile(baseDirs{bi},subdirs{bi}{di},'filtered',sprintf(bare,coli-1,chanid-1));
        im = imread(fname);
        lim = stitchedlim(im,tols(ci,:));
        lim(2) = lim(1) + scale*(lim(2) - lim(1));
        lims(bi,:) = lim;
    end

    ims = cell(1,N);
    for pidx = 1:N
        coli = piecols(pidx);
        bi = piedirs(pidx);
        cii = find(strcmp(channelLabels{bi},chan));
        di = channelDirs{bi}(cii);
        chanid = channelIds{bi}(cii);
        %load image
        fname = fullfile(baseDirs{bi},subdirs{bi}{di},'filtered',sprintf(bare,coli-1,chanid-1));
        im = imread(fname);
        %adjust contrast, center, & crop
        im = im2double(imadjust(im,lims(bi,:)));
        %im = mat2gray(im);
        if pidx == 3
            im = imresize(im, [2765,3156.9]);
        end

        % centered crop indices
        if pidx == 3
            CM = round(pos{bi}(coli).center*(2765/1277));% + margin;
        else
            CM = round(pos{bi}(coli).center);% + margin;
        end

        if pidx == 2
            CM = CM+[0,50];
        end
        if pidx == 3
            CM = CM+[0,50];
        end
        if pidx == 4
            CM = CM+[0,50];
        end
        %CM = round(pos{bi}(coli).center);% + margin;
        yrange = (CM(2)-Rmax):(CM(2)+Rmax);
        xrange = (CM(1)-Rmax):(CM(1)+Rmax);

        im = im(yrange, xrange);
        I = ones(cropsize,cropsize);
        I(padsize+1:padsize+1+2*Rmax,padsize+1:padsize+1+2*Rmax) = im;
        I(mask{pidx}) = 0;
        I(disk) = 1;
        I(lines) = 1;

        ims{pidx} = I;
    end

    combined = sum(cat(3,ims{:}),3);
    clf
    imshow(combined)
    title(chan)
    drawnow

    savename = fullfile(savedir,[pref,'_',chan,'.png']);
    imwrite(combined,savename)
end












