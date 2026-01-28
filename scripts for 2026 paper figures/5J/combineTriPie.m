clc
clear
file = 'combinedRadialSOX17.xlsx';
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
sheetList = sheetnames(file);
for k = 1:numel(sheetList)
    tbl = readtable(file,'Sheet',sheetList{k});
    res{k}.r = tbl.r; % Distance
    res{k}.nuc_profile = table2array(tbl(:,startsWith(tbl.Properties.VariableNames,'nuc_profile')));
    res{k}.nuc_profile = res{k}.nuc_profile(:,3); % both SOX17 and TBX6 are in channel 3
    res{k}.nuc_profile_colstd = table2array(tbl(:,startsWith(tbl.Properties.VariableNames,'nuc_profile_colstd')));
    res{k}.nuc_profile_colstd = res{k}.nuc_profile_colstd(:,3);
end

conditionindices = [1,3,4];
res = res(conditionindices);
colors = lines(numel(conditionindices));
lw = 3;
fs = 28;

%markerChannels = options.nucChannels+1;
for channel = 1
    figure('Position',[0 0 600 280])
    legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = 1:numel(conditionindices)
        i = i+1;

        plot(res{condi}.r, res{condi}.nuc_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));

    end
    i = 0;
    for condi = 1:numel(conditionindices)
        i = i+1;

        plerr = res{condi}.nuc_profile(:,channel) + res{condi}.nuc_profile_colstd(:,channel);
        minerr = res{condi}.nuc_profile(:,channel) - res{condi}.nuc_profile_colstd(:,channel);
        good = ~isnan(plerr);
        xvals = res{condi}.r(good);
        y_upper = plerr(good);
        y_lower = minerr(good);

        % Ensure row vectors
        xvals = xvals(:)';
        y_upper = y_upper(:)';
        y_lower = y_lower(:)';

        XX = [xvals, fliplr(xvals)];
        YY = [y_upper, fliplr(y_lower)];

        fill(XX, YY, colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    hold off

    xlim([0 350]);
    ylim([0 1.25]);

    ylabel('intensity');

    if numel(conditionindices) == 4
    legendstr = {'ctrl shRNA','ctrl shRNA','FGF17 shRNA','FGF4 shRNA'};
    end
    if numel(conditionindices) == 3
    legendstr = {'ctrl shRNA','FGF17 shRNA','FGF4 shRNA'};
    end
    legend(p,legendstr,'FontSize',10,'Location','Northeast')

    set(gca,'FontSize', fs)
    xlabel('edge distance ( um )','FontSize',40)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', lw);
    %axis square;
    saveas(gcf, fullfile(scriptPath, ['radialintensity_combine_tri pie SOX17.png']));
end