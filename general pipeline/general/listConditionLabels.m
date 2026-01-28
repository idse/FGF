function condlabels = listConditionLabels(meta)

npos = meta.nPositions;
conditionStartPos = meta.conditionStartPos;

if npos ~= sum(meta.posPerCondition)
    error('posPerCondition does not match nPositions')
end

condlabels = cell(1,npos);
for ii = 1:npos
    condi = find(ii >= conditionStartPos,1,'last'); % condition index
    condlabels{ii} = meta.conditions{condi};
end


end