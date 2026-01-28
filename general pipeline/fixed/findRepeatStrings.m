function [repeatStrings,poslabels] = findRepeatStrings(barenames)

s = strsplit(barenames{1},'_');
repeatStrings = cell(size(s));
repeatMask = false(size(s));
for ii = 1:length(s)
    pattern = ['_',s{ii},'_'];
    if ii == 1
        pattern = pattern(2:end);
    end
    
    if ii == length(s)
        pattern = pattern(1:end-1);
    end
    
    if all(contains(barenames,pattern))
        repeatStrings{ii} = pattern;
        repeatMask(ii) = true;
    end
end

cc = bwconncomp(repeatMask);
idxs = cc.PixelIdxList;
newRepeats = cell(size(idxs));
for ii = 1:length(idxs)
    t = repeatStrings(idxs{ii});
    newRepeats{ii} = strrep(strcat(t{:}),'__','_');
end
repeatStrings = newRepeats;

poslabels = cell(size(barenames));
for ii = 1:length(barenames)
    name = barenames{ii};
    oldname = name;
    for jj = 1:length(repeatStrings)
        s = strsplit(oldname,repeatStrings{jj});
        if isempty(s{1}) || isempty(s{end})
            name = strrep(name,repeatStrings{jj},'');
        else
            name = strrep(name,repeatStrings{jj},'_');
        end
    end
    poslabels{ii} = name;
end

end