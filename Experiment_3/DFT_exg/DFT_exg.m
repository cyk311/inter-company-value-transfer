%% exg_plot_dft.m
% 任务：
% 1) 一张大图（上下两幅子图）
% 2) 上图：a 为折线；a_hat 为平滑曲线
% 3) 下图：对 a 做 DFT，绘制 Magnitude-Frequency，并标注最大强度频率（显示取整）

clear; clc; close all;

%% 读取数据
file = 'exg.xlsx';
T = readtable(file);

% 基本列检查（若列名不一致，请在这里改）
year  = T.year;
a     = T.a;
a_hat = T.a_hat;

% 排序 + 去缺失（保证 FFT 输入干净且等间隔的年序列）
TT = table(year, a, a_hat);
TT = sortrows(TT, "year");
TT = rmmissing(TT);

year  = TT.year;
a     = TT.a;
a_hat = TT.a_hat;

%% 1) 生成图像框：一张大图
fig = figure('Color','w', 'Units','pixels', 'Position',[100 100 1400 900], ...
             'Name','a vs a\_hat and DFT(a)', 'NumberTitle','off');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

%% 2) 子图一（上侧）：a 折线 + a_hat 平滑曲线
ax1 = nexttile;
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

% a：折线
plot(ax1, year, a, '-', 'LineWidth', 1.4, 'DisplayName','$a$');

% a_hat：平滑曲线（用 pchip 插值到更密集的年点，得到“平滑曲线”观感）
yearFine = linspace(min(year), max(year), max(200, numel(year)*10));
aHatFine = interp1(year, a_hat, yearFine, 'pchip');
plot(ax1, yearFine, aHatFine, 'LineWidth', 2.2, 'DisplayName','$\hat{a}$ (smooth)');
lgd = legend(ax1, 'show');
lgd.Interpreter = 'latex';
lgd.FontSize = 13;

xlabel(ax1, 'Year');
ylabel(ax1, 'Value');
title(ax1, 'Exg. $a$ and $\hat{a}$ series when $i=16$ \& $j=5$', 'Interpreter','latex','FontSize',13);
legend(ax1, 'Location','best');

%% 3) 子图二（下侧）：对 a 做 DFT，绘制 Magnitude-Frequency，并标注最大频率（取整显示）
ax2 = nexttile;
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

x = a;

% 常规做法：先去均值，避免 0 频（DC）分量压制其他频率
x = x - mean(x, 'omitnan');

N = numel(x);
Fs = 1; % 年度数据：采样频率 = 1 次/年

Y = fft(x);
P2 = abs(Y) / N;                         % 双边幅度
P1 = P2(1:floor(N/2)+1);                 % 单边幅度
if numel(P1) > 2
    P1(2:end-1) = 2*P1(2:end-1);         % 单边谱能量补偿（除 DC 与 Nyquist）
end
f = (0:floor(N/2))*(Fs/N);               % 频率轴（cycles/year）

plot(ax2, f, P1, 'LineWidth', 1.4);
xlabel(ax2, 'Frequency (cycles/year)');
ylabel(ax2, 'Magnitude');
title(ax2, 'DFT Magnitude Spectrum of $a$ (demeaned)', 'Interpreter','latex','FontSize',13);

% 找最大强度频率（忽略 f=0 的 DC 项）
if numel(P1) >= 2
    [~, idxRel] = max(P1(2:end));
    idxPeak = idxRel + 1;
else
    idxPeak = 1;
end

fPeak   = f(idxPeak);
magPeak = P1(idxPeak);
fPeakInt = round(fPeak);  % “取整”用于显示

% 标注峰值
plot(ax2, fPeak, magPeak, 'o', 'MarkerSize', 8, 'LineWidth', 1.6);
xline(ax2, fPeak, '--', 'LineWidth', 1.2);

txt = sprintf('Peak f = %.4f (rounded: %d)', fPeak, fPeakInt);
text(ax2, fPeak, magPeak, ['  ' txt], 'VerticalAlignment','bottom', 'FontSize', 11);

% 年度数据 Nyquist 频率为 0.5 cycles/year
xlim(ax2, [0, 0.5]);

%% （可选）导出图片
% exportgraphics(fig, 'exg_a_ahat_dft.png', 'Resolution', 300);
