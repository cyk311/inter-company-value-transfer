%% count_plot_Nt.m
% Requirements:
%   value_transfer_long.csv contains columns: year, ue_ij

clear; clc;

file = 'value_transfer_long.csv';
assert(isfile(file), 'Cannot find %s in current folder.', file);

T = readtable(file);

% --- Robustly locate columns ---
vars = lower(string(T.Properties.VariableNames));

iy = find(vars=="year", 1);
iu = find(vars=="ue_ij" | vars=="ueij" | vars=="ue", 1);

assert(~isempty(iy), 'Column "year" not found in %s.', file);
assert(~isempty(iu), 'Column "ue_ij" (or ueij/ue) not found in %s.', file);

year = double(T{:,iy});
ue   = double(T{:,iu});

% Drop rows with missing year
maskY = ~isnan(year);
year = year(maskY);
ue   = ue(maskY);

% --- Compute n^t, m^t, N^t by year ---
years = unique(year);
years = sort(years);

n = zeros(numel(years),1);
m = zeros(numel(years),1);
N = nan(numel(years),1);

for k = 1:numel(years)
    yk = years(k);
    idx = (year==yk);

    ue_k = ue(idx);
    valid = ~isnan(ue_k);

    m(k) = sum(valid);
    n(k) = sum(valid & (ue_k > 1));

    if m(k) > 0
        N(k) = n(k) / m(k);
    end
end

% Optional: output yearly summary table
summaryTbl = table(years, n, m, N, 'VariableNames', {'year','n_t','m_t','N_t'});
writetable(summaryTbl, 'Nt_summary.csv');

% --- Plot N^t time series with reference line at 0.5 ---
figure('Color','w');
plot(years, N, 'LineWidth', 1.8); hold on;
yline(0.5, 'r-', 'LineWidth', 1.6);  % red reference line

grid on;
xlim([min(years) max(years)]);
ylim([0 1]); % optional; comment out if you prefer auto
xlabel('Year');
ylabel('N^t = n^t / m^t');
title('Share of UE entries greater than 1 by year');

% If you want to save the figure:
% exportgraphics(gcf, 'Nt_timeseries.png', 'Resolution', 300);

disp('Done.');
disp('Wrote: Nt_summary.csv');
