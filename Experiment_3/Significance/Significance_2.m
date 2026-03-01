%% plot_ue_significance_heatmap_gradient.m
clear; clc;

%% ===== Settings =====
xlsxFile   = 'ue_significance_matrix.xlsx';
sheetClass = 'class';
sheetP     = 'p_value';

n     = 18;
alpha = 0.05;
rangeMat = 'B2:S19';

p_floor = 1e-6;
use_log_scale = true;

show_p_in_nonsig_only = true;
p_fmt = '%.3f';
p_small_as = '<0.001';

gridColor = [0.85 0.85 0.85];
gridLW = 0.8;

outPng = 'ue_significance_heatmap_gradient.png';
outPdf = 'ue_significance_heatmap_gradient.pdf';
export_pdf = true;
% =====================

assert(isfile(xlsxFile), 'Cannot find %s', xlsxFile);

classMat = readmatrix(xlsxFile, 'Sheet', sheetClass, 'Range', rangeMat);
pMat     = readmatrix(xlsxFile, 'Sheet', sheetP,     'Range', rangeMat);

if ~isequal(size(classMat), [n n]) || ~isequal(size(pMat), [n n])
    error('Matrix sizes are not %dx%d. Check sheets/range.', n, n);
end

%% Build intensity matrix C in [-1,1]
C = zeros(n,n);

pEff = pMat;
pEff(~isfinite(pEff)) = 1;
pEff = max(pEff, p_floor);

if use_log_scale
    Lalpha = -log10(alpha);
    Lfloor = -log10(p_floor);
    inten  = (-log10(pEff) - Lalpha) ./ (Lfloor - Lalpha);
else
    inten = 1 - (pEff./alpha);
end
inten = max(0, min(1, inten));

maskIn  = (classMat== 1);
maskOut = (classMat==-1);

C(maskIn)  =  inten(maskIn);
C(maskOut) = -inten(maskOut);
C(~isfinite(classMat)) = 0;

%% Colormap: continuous blue-white-red
m = 256; half = floor(m/2);
blue  = [0 0 1]; white = [1 1 1]; red = [1 0 0];
cmap1 = [linspace(blue(1),  white(1), half)', ...
         linspace(blue(2),  white(2), half)', ...
         linspace(blue(3),  white(3), half)'];
cmap2 = [linspace(white(1), red(1),   m-half)', ...
         linspace(white(2), red(2),   m-half)', ...
         linspace(white(3), red(3),   m-half)'];
cmap  = [cmap1; cmap2];

%% Plot (FIXED alignment)
figure('Color','w','Position',[120 80 980 760]);
ax = axes();

% 关键修正：用 1:n 作为像素中心，避免 0.5 格偏移
imagesc(ax, 1:n, 1:n, C);

axis(ax, 'image');
set(ax, 'YDir','reverse');                 % 让 i=1 在最上方
xlim(ax, [0.5 n+0.5]);
ylim(ax, [0.5 n+0.5]);

set(ax, 'XTick', 1:n, 'YTick', 1:n, 'TickDir','out');
ax.TickLength = [0 0];

xlabel(ax, 'Downstream sector j');
ylabel(ax, 'Upstream sector i');
title(ax, sprintf('UE direction significance (alpha=%.2f): red=inflow, blue=outflow; shade by p', alpha));

colormap(ax, cmap);
caxis(ax, [-1 1]);

% Gridlines on cell borders (half-integers)
hold(ax, 'on');
for k = 0.5:1:(n+0.5)
    line(ax, [0.5 n+0.5], [k k], 'Color', gridColor, 'LineWidth', gridLW);
    line(ax, [k k], [0.5 n+0.5], 'Color', gridColor, 'LineWidth', gridLW);
end

% Annotate p-values in non-significant cells
for ii = 1:n
    for jj = 1:n
        if ~isfinite(classMat(ii,jj)) || ~isfinite(pMat(ii,jj)), continue; end
        isNS = (classMat(ii,jj)==0);
        if show_p_in_nonsig_only && ~isNS, continue; end

        p = pMat(ii,jj);
        if p < 0.001
            txt = p_small_as;
        else
            txt = sprintf(p_fmt, p);
        end

        text(ax, jj, ii, txt, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',8, ...
            'Color',[0 0 0]);
    end
end
hold(ax, 'off');

cb = colorbar(ax);
cb.Ticks = [-1 -0.5 0 0.5 1];
cb.TickLabels = {'Strong outflow','Weak outflow','Not significant','Weak inflow','Strong inflow'};

exportgraphics(gcf, outPng, 'Resolution', 300);
fprintf('Saved: %s\n', outPng);

if export_pdf
    exportgraphics(gcf, outPdf, 'ContentType','vector');
    fprintf('Saved: %s\n', outPdf);
end
