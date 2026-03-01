%% plot_ue_3x3_logratio.m
% 3x3 tiles for selected years.
% Each tile: left = log(UE_i) strip, right = log(UE_ij) matrix
% Diverging colormap centered at 0: blue (neg) -> white (0) -> red (pos)
% Unified colorbar on the right.

clear; clc;

%% ===== User settings =====
inFile = 'value_transfer_long.csv';
yearsPlot = [1986 1990 1994 1998 2002 2006 2010 2014 2018];
n = 18;

% Layout controls
rightMargin = 0.095;
tilePad     = 0.012;
vecWidth    = 0.05;      % UE_i strip width
gap         = 0.03;
titleBand   = 0.085;
fontSize    = 9;

% Extra gaps between tiles (continuous control)
tileGapX = 0.006;   % horizontal half-gap per tile
tileGapY = 0.008;   % vertical half-gap per tile


% Log-ratio settings
logBase   = "ln";        % "ln" or "log10"
epsUE     = 1e-12;       % floor UE to avoid log(0)
clipQuant = 0.035;       % robust limits in transformed space
nColors   = 257;         % odd -> exact white at center
boostK    = 6;           % <<< moved here (was inside loop)

saveFig = false;
outPng  = 'UE_3x3_logratio.png';
dpi     = 300;
%% =========================

assert(isfile(inFile), 'Cannot find %s', inFile);
T = readtable(inFile);

vars = lower(string(T.Properties.VariableNames));
req = ["year","i","j","ue_ij"];
for k = 1:numel(req)
    assert(any(vars==req(k)), 'CSV missing required column: %s', req(k));
end

year = T{:, vars=="year"};
i    = T{:, vars=="i"};
j    = T{:, vars=="j"};
ueij = T{:, vars=="ue_ij"};

% Optional UE_i column
hasUEi = any(vars=="ue_i") || any(vars=="uei") || any(vars=="ue_i_vector");
if hasUEi
    if any(vars=="ue_i");       uei_col = T{:, vars=="ue_i"};
    elseif any(vars=="uei");    uei_col = T{:, vars=="uei"};
    else;                       uei_col = T{:, vars=="ue_i_vector"};
    end
else
    uei_col = nan(size(ueij));
end

maskY = ismember(year, yearsPlot);
assert(any(maskY), 'No rows match requested years in %s.', inFile);

year = year(maskY); i = i(maskY); j = j(maskY); ueij = ueij(maskY); uei_col = uei_col(maskY);
assert(all(i>=1 & i<=n & j>=1 & j<=n), 'i/j out of range 1..%d', n);

%% Helper: log transform with floor
logfun = @(x) local_log_ratio(x, epsUE, logBase);

%% Build per-year UE_i and UE_ij (then transform), and compute global color limits
allVals = [];
yearData = struct();

for t = 1:numel(yearsPlot)
    yr = yearsPlot(t);
    idx = (year==yr);

    UEij = nan(n,n);
    UEi  = nan(n,1);

    if any(idx)
        ii = i(idx); jj = j(idx); ue = ueij(idx);
        UEij(sub2ind([n n], ii, jj)) = ue;

        if hasUEi
            uei_this = uei_col(idx);
            for s = 1:n
                UEi(s) = mean(uei_this(ii==s), 'omitnan');
            end
        else
            for s = 1:n
                UEi(s) = mean(ue(ii==s), 'omitnan'); % fallback
            end
        end
    end

    % log-ratio then non-linear boost
    LUEi  = asinh(boostK * logfun(UEi));
    LUEij = asinh(boostK * logfun(UEij));

    yearData(t).yr    = yr;
    yearData(t).LUEi  = LUEi;
    yearData(t).LUEij = LUEij;

    allVals = [allVals; LUEi(:); LUEij(:)]; %#ok<AGROW>
end

allVals = allVals(isfinite(allVals));
assert(~isempty(allVals), 'No finite transformed values found for selected years.');

% Robust limits, force symmetry around 0
allVals = sort(allVals);
qlo = quantile_no_toolbox(allVals, clipQuant);
qhi = quantile_no_toolbox(allVals, 1-clipQuant);

d = max(abs(qlo), abs(qhi));
if ~isfinite(d) || d==0, d = 1; end
clim = [-d, d];

%% Figure + tiled layout (3x3)
fig = figure('Color','w');
fig.Position(3) = 1400;
fig.Position(4) = 900;

tl = tiledlayout(fig, 3, 3, 'TileSpacing','tight', 'Padding','compact');
tl.OuterPosition = [0 0 1-rightMargin 1];

cmap = diverging_bwr(nColors);
colormap(fig, cmap);

for k = 1:numel(yearsPlot)

    % --- IMPORTANT FIX ---
    % Use a "tile anchor axes" to get a consistent tile rectangle.
    % Do NOT use ax.Position directly without normalizing insets.
    axTmp = nexttile(tl, k);
    set(axTmp, 'Units','normalized', 'Visible','off', ...
        'XTick',[], 'YTick',[], 'XTickLabel',[], 'YTickLabel',[], ...
        'Box','off', 'LooseInset',[0 0 0 0]);

    drawnow limitrate;
    tilePos = axTmp.OuterPosition;   % <<< changed from axTmp.Position
    tilePos = tilePos + [tileGapX, tileGapY, -2*tileGapX, -2*tileGapY];

    % Safety (avoid negative widths/heights if you set too large gaps)
    tilePos(3) = max(tilePos(3), 0.01);
    tilePos(4) = max(tilePos(4), 0.01);
    % ----------------------

    p = uipanel(fig, 'Units','normalized', 'Position', tilePos, 'BorderType','none');

    % Title band
    uicontrol('Parent', p, 'Style','text', 'Units','normalized', ...
        'Position', [0, 1-titleBand, 1, titleBand], ...
        'String', sprintf('%d', yearData(k).yr), ...
        'FontWeight','bold', 'FontSize', 11, ...
        'BackgroundColor','w', 'HorizontalAlignment','center');

    x0 = tilePad; y0 = tilePad;
    w0 = 1 - 2*tilePad;
    h0 = 1 - 2*tilePad - titleBand;

    wVec = w0 * vecWidth;
    wMat = w0 - wVec - gap;

    axVec = axes('Parent', p, 'Units','normalized', ...
        'Position', [x0, y0, wVec, h0]);
    axMat = axes('Parent', p, 'Units','normalized', ...
        'Position', [x0 + wVec + gap, y0, wMat, h0]);

    LUEi  = yearData(k).LUEi;
    LUEij = yearData(k).LUEij;

    imagesc(axVec, LUEi);
    set(axVec, 'YDir','normal', 'XTick',[], 'YTick',[], ...
        'FontSize', fontSize, 'TickLength',[0 0]);
    ylim(axVec, [0.5 n+0.5]);
    box(axVec, 'on');
    caxis(axVec, clim);

    imagesc(axMat, LUEij);
    set(axMat, 'YDir','normal', 'XTick',[], 'YTick',[], ...
        'FontSize', fontSize, 'TickLength',[0 0]);
    axis(axMat, 'tight');
    box(axMat, 'on');
    caxis(axMat, clim);
end

%% Unified colorbar (dummy axis)
axCB = axes(fig, 'Units','normalized', ...
    'Position', [1-rightMargin+0.02, 0.12, 0.001, 0.76], ...
    'Visible','off');
colormap(axCB, cmap);
caxis(axCB, clim);

cb = colorbar(axCB);
cb.Units = 'normalized';
cb.Position = [1-rightMargin+0.045, 0.12, 0.02, 0.76];
cb.FontSize = fontSize;

% --- label FIX: ln vs log10 was reversed in your script ---
if logBase=="log10"
    cb.Label.String = sprintf('asinh(%g·log_{10}(UE))', boostK);
else
    cb.Label.String = sprintf('asinh(%g·ln(UE))', boostK);
end
cb.Label.FontSize = 11;

cb.Ticks = unique([clim(1), 0, clim(2)]);
cb.TickLabels = compose('%.2f', cb.Ticks);

drawnow;

if saveFig
    exportgraphics(fig, outPng, 'Resolution', dpi);
    fprintf('Saved: %s\n', outPng);
end

%% ---- local helpers ----
function y = local_log_ratio(x, epsUE, logBase)
    x = double(x);
    y = nan(size(x));
    m = isfinite(x);
    if any(m(:))
        x2 = x;
        x2(m) = max(x2(m), epsUE);
        if logBase=="log10"
            y(m) = log10(x2(m));
        else
            y(m) = log(x2(m));
        end
    end
end

function q = quantile_no_toolbox(sortedVec, p)
    n = numel(sortedVec);
    if n==0 || ~isfinite(p), q = NaN; return; end
    p = max(0, min(1, p));
    idx = 1 + (n-1)*p;
    i0 = floor(idx); i1 = ceil(idx);
    i0 = max(1, i0); i1 = min(n, i1);
    if i0==i1
        q = sortedVec(i0);
    else
        w = idx - i0;
        q = (1-w)*sortedVec(i0) + w*sortedVec(i1);
    end
end

function cmap = diverging_bwr(n)
    n = max(3, round(n));
    if mod(n,2)==0, n = n + 1; end
    m = (n+1)/2;

    blue  = [0.02, 0.10, 0.75];
    white = [1, 1, 1];
    red   = [0.75, 0.02, 0.05];

    a = linspace(0,1,m)';  left  = (1-a).*blue  + a.*white;
    b = linspace(0,1,m)';  right = (1-b).*white + b.*red;

    cmap = [left(1:end-1,:); white; right(2:end,:)];
    cmap(m,:) = white;
end
