function res = makeSegmentationCrossSection(img,LM,meta,yinds,xinds)

%img = nuclear channel image stack
%L = 3D label matrix
%meta = experiment metadata
%xinds = x indices
%yinds = y indices

% if numel(xinds) > 1 && numel(yinds) > 1
%     error('either x or y range must be only a single index')
% end

if numel(xinds) < numel(yinds)
    disp('x section')    
    res = meta.yres;
    im = squeeze(max(img(yinds,xinds,:),[],2))';
    L = squeeze(max(LM(yinds,xinds,:),[],2))';
else
    disp('y section')
    res = meta.xres;
    im = squeeze(max(img(yinds,xinds,:),[],1))';
    L = squeeze(max(LM(yinds,xinds,:),[],1))';
end

im = imadjust(im);
segoverlay = visualize_nuclei_v2(L,im);

ratio = meta.zres/res;
newsize = round([size(im,1)*ratio size(im,2)]);

im = flipud(imresize(im,newsize));
segoverlay = flipud(imresize(segoverlay,newsize));

figure
subplot_tight(2,1,1)
imshow(im)
title('Image')
%cleanSubplot

subplot_tight(2,1,2)
imshow(segoverlay)
title('With Labels')
%cleanSubplot

res = struct('xsection', im, 'segoverlay', segoverlay);

end