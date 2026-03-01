%% rho_matrix_1981_2018_plot.m
clear; clc; close all;

file = 'rho_matrix_1981_2018.xlsx';
thr  = 1;

% ---------- Read data ----------
T = readtable(file, 'VariableNamingRule','preserve');
vnames = string(T.Properties.VariableNames);

yearIdx = find(lower(vnames)=="year", 1);
if ~isempty(yearIdx)
    year = T{:, yearIdx};
    dataCols = setdiff(1:width(T), yearIdx);
else
    firstCol = T{:,1};
    if isnumeric(firstCol) && all(firstCol>=1900 & firstCol<=2100)
        year = firstCol;
        dataCols = 2:width(T);
    else
        year = (1981:2018).';
        dataCols = 1:width(T);
    end
end

Y = T{:, dataCols};                 % nYear x 18
if size(Y,2) ~= 18
    warning('Detected %d series columns (expected 18). Please check the Excel layout.', size(Y,2));
end

[year, order] = sort(year);
Y = Y(order, :);

% ---------- Sector names (01-18) ----------
sectorNames = { ...
'农林牧渔产品和服务'
'采掘产品'
'食品和烟草'
'纺织、服装、鞋帽及皮革羽绒制品'
'木材加工、家具、造纸印刷和文教体育用品'
'石油、炼焦、核燃料加工品和化学产品'
'非金属矿物制品'
'金属冶炼、加工及制品'
'机械设备、交通运输设备、电子电气及其他设备'
'其他制造产品及修理服务'
'电力、热力、燃气及水的生产和供应'
'建筑'
'批发和零售'
'交通运输、仓储和邮政'
'信息传输、软件和信息技术服务'
'金融和房地产'
'研究和试验发展'
'其他服务' };

% ---------- Plot ----------
fig = figure('Color','w');
n = size(Y,2);

tlo = tiledlayout(fig, n, 1, 'TileSpacing','none', 'Padding','compact');

% 给左侧名称留出空间（关键：避免被 figure 边界裁切）
tlo.OuterPosition = [0.25 0.02 0.78 0.96];   % 左侧留空间，图表贴靠右侧

ax = gobjects(n,1);

for k = 1:n
    ax(k) = nexttile(tlo);
    y = Y(:,k);

    hold(ax(k),'on');
    yline(ax(k), thr, '--', 'LineWidth', 0.8);
    plotColoredByThreshold(ax(k), year, y, thr);

    grid(ax(k),'on');
    ax(k).Box = 'on';

    % 不要 y 轴 tick / tick label
    set(ax(k), 'YTick', [], 'YTickLabel', []);

    % 左侧写部门名称（不作为 y-label；放在轴外侧）
    txt = sectorNames{k};
    text(ax(k), -0.02, 0.5, txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontSize', 9, 'FontName','Microsoft YaHei', ... % 可改 SimHei 等
        'Clipping','off', 'Interpreter','none');

    % y-limits with margin
    yy = y(~isnan(y));
    if ~isempty(yy)
        yMin = min(yy); yMax = max(yy);
        if yMin == yMax
            pad = max(0.05*abs(yMin), 0.1);
        else
            pad = 0.08*(yMax - yMin);
        end
        ylim(ax(k), [yMin-pad, yMax+pad]);
    end

    if k < n
        ax(k).XTickLabel = [];
    else
        xlabel(ax(k), 'Year');
    end

    hold(ax(k),'off');
end

linkaxes(ax, 'x');
xlim(ax(1), [min(year), max(year)]);

% exportgraphics(fig, 'rho_18series_piecewise_colored.png', 'Resolution', 300);

%% ---------- Local function ----------
function plotColoredByThreshold(ax, x, y, thr)
    blue = 'b';
    red  = 'r';
    lw   = 1.2;

    x = x(:); y = y(:);
    n = numel(x);

    for i = 1:n-1
        x1 = x(i);   x2 = x(i+1);
        y1 = y(i);   y2 = y(i+1);

        if any(isnan([x1 x2 y1 y2])), continue; end

        if y1 == thr && y2 == thr
            line(ax, [x1 x2], [y1 y2], 'Color', blue, 'LineWidth', lw);
            continue;
        end

        if (y1 <= thr && y2 <= thr)
            line(ax, [x1 x2], [y1 y2], 'Color', blue, 'LineWidth', lw);
        elseif (y1 >= thr && y2 >= thr)
            line(ax, [x1 x2], [y1 y2], 'Color', red, 'LineWidth', lw);
        else
            t  = (thr - y1) / (y2 - y1);
            xC = x1 + t*(x2 - x1);

            if y1 < thr
                line(ax, [x1 xC], [y1 thr], 'Color', blue, 'LineWidth', lw);
                line(ax, [xC x2], [thr y2], 'Color', red,  'LineWidth', lw);
            else
                line(ax, [x1 xC], [y1 thr], 'Color', red,  'LineWidth', lw);
                line(ax, [xC x2], [thr y2], 'Color', blue, 'LineWidth', lw);
            end
        end
    end
end
