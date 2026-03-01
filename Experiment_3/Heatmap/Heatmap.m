%% plot_ue_heatmaps_by_year.m
% Generate 38 yearly heatmaps (18x18) for ue_ij from value_transfer_long.csv

clear; clc;

% -------- User settings --------
inFile   = 'value_transfer_long.csv';
outDir   = 'ue_heatmaps';
saveFIG  = false;   % set true if you also want .fig files
useGlobalCaxis = true;  % keep the same color scale across years
useRobustCaxis = true;  % use percentile clipping to reduce outlier domination
robustLo = 1;          % lower percentile
robustHi = 99;         % upper percentile
% -------------------------------

assert(isfile(inFile), 'Cannot find %s', inFile);
if ~exist(outDir, 'dir'), mkdir(outDir); end

T = readtable(inFile);

% Basic checks
needVars = {'year','i','j','ue_ij'};
assert(all(ismember(needVars, T.Properties.VariableNames)), ...
    'Input must contain columns: year, i, j, ue_ij');

T.year = double(T.year);
T.i    = double(T.i);
T.j    = double(T.j);

years = unique(T.year);
years = sort(years(:));
n = 18;

% Optional axis labels from sector_map.csv (if exists)
xLabels = string(1:n);
yLabels = string(1:n);
if isfile('sector_map.csv')
    sm = readtable('sector_map.csv');
    if all(ismember({'sector_id','sector_name'}, sm.Properties.VariableNames))
        sm = sortrows(sm,'sector_id');
        if height(sm) >= n
            names = string(sm.sector_name(1:n));
            xLabels = names;
            yLabels = names;
        end
    end
end

% Compute global color limits (optional)
if useGlobalCaxis
    ueAll = T.ue_ij;
    ueAll = ueAll(isfinite(ueAll));
    if isempty(ueAll)
        error('No finite ue_ij values found.');
    end
    if useRobustCaxis
        cmin = prctile(ueAll, robustLo);
        cmax = prctile(ueAll, robustHi);
    else
        cmin = min(ueAll);
        cmax = max(ueAll);
    end
    if cmin == cmax
        cmin = cmin - 1;
        cmax = cmax + 1;
    end
end

fprintf('Generating %d heatmaps into folder: %s\n', numel(years), outDir);

for yy = 1:numel(years)
    yr = years(yy);
    S = T(T.year==yr, :);

    % Build 18x18 matrix with NaN default
    UE = nan(n,n);
    idx = ~isnan(S.i) & ~isnan(S.j) & S.i>=1 & S.i<=n & S.j>=1 & S.j<=n;
    ii = S.i(idx);
    jj = S.j(idx);
    vv = S.ue_ij(idx);

    % If duplicates exist, take mean
    UE = accumarray([ii, jj], vv, [n n], @mean, NaN);

    % Plot
    f = figure('Visible','off', 'Color','w');
    ax = axes(f);

    h = imagesc(ax, UE);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');

    % Make NaNs transparent (so background shows through)
    set(h, 'AlphaData', isfinite(UE));

    colormap(ax, parula);
    cb = colorbar(ax);
    cb.Label.String = 'UE_{ij} = p_{ij} / \lambda_i';

    title(ax, sprintf('UE heatmap (ue_{ij}) - %d', yr), 'Interpreter','none');
    xlabel(ax, 'Downstream sector j');
    ylabel(ax, 'Upstream sector i');

    xticks(ax, 1:n); yticks(ax, 1:n);
    xticklabels(ax, xLabels); yticklabels(ax, yLabels);
    xtickangle(ax, 45);

    if useGlobalCaxis
        caxis(ax, [cmin cmax]);
    end

    % Save
    outPng = fullfile(outDir, sprintf('UE_%d.png', yr));
    exportgraphics(f, outPng, 'Resolution', 200);

    if saveFIG
        savefig(f, fullfile(outDir, sprintf('UE_%d.fig', yr)));
    end

    close(f);

    if mod(yy,10)==0 || yy==numel(years)
        fprintf('  %d / %d done\n', yy, numel(years));
    end
end

fprintf('All done.\n');
