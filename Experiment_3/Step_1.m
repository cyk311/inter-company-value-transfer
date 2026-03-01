%% 第一部分：将投入产出表转化为面板数据格式
% 读取IO表 (1981-2018), 计算技术系数 a_ij^t = x_ij^t / x_j,
% 导出为面板数据格式: year,i,j,a

clear; clc;

%% ===== USER SETTINGS =====
ioFile   = 'IO(1981-2018).xlsx';     % your workbook
outDir   = '.';                      % output folder
rangeNum = 'D9:AI33';                % numeric block (25x32) per your definition
rangeRowNames = 'B9:B26';            % row names (18 industries + totals etc.)
rangeColNames = 'D5:U5';             % column names (18 industries)
yearsWanted = (1981:2018)';          % 38 years
n = 18;
% ==========================

assert(isfile(ioFile), 'Cannot find %s. Put it in the working directory or change ioFile.', ioFile);

% Detect available sheets and intersect with requested years
allSheets = sheetnames(ioFile);
yrAvail = [];
for k = 1:numel(allSheets)
    yk = str2double(allSheets{k});
    if ~isnan(yk)
        yrAvail(end+1,1) = yk; %#ok<SAGROW>
    end
end
years = intersect(yearsWanted, unique(yrAvail));
assert(~isempty(years), 'No year-named sheets found that match 1981-2018.');

T = numel(years);
Nobs = n*n*T;

yearVec = zeros(Nobs,1,'int32');
iVec    = zeros(Nobs,1,'int16');
jVec    = zeros(Nobs,1,'int16');
aVec    = nan(Nobs,1);
zVec    = nan(Nobs,1);
xjVec   = nan(Nobs,1);

% Read sector names (take first available year sheet to extract labels)
refSheet = num2str(years(1));
rowNames = readcell(ioFile,'Sheet',refSheet,'Range',rangeRowNames);
colNames = readcell(ioFile,'Sheet',refSheet,'Range',rangeColNames);
rowNames = rowNames(1:n);
colNames = colNames(1:n);

% Export sector map (optional but useful)
sector_id = (1:n)';
sector_name = string(rowNames);
writetable(table(sector_id, sector_name), fullfile(outDir,'sector_map.csv'));

% Precompute index grids (column-major friendly: j outer, i inner)
[I,J] = ndgrid(1:n,1:n);
I = I(:); J = J(:); % i changes fastest within each j? Actually ndgrid gives i varying fastest within column-major fill.

ptr = 0;

for tt = 1:T
    yr = years(tt);
    sh = num2str(yr);

    M = readmatrix(ioFile,'Sheet',sh,'Range',rangeNum);

    % Extract intermediate transactions Z (18x18)
    Z = M(1:n, 1:n);

    % Gross output x (use row total "总产出" column=32, rows 1..18)
    x = M(1:n, 32);

    if any(~isfinite(x)) || any(x==0)
        warning('Year %d: x has NaN/Inf/0. A_ij will be NaN where x_j invalid.', yr);
    end

    % Direct requirement matrix A = Z * diag(1./x)
    % (divide each column j by x_j)
    A = Z ./ (x');  % implicit expansion (R2016b+). For older MATLAB, use bsxfun(@rdivide,Z,x')

    vecA = A(:);
    vecZ = Z(:);
    vecXj = repmat(x(:), n, 1); % stacks x for each i within each j

    idx = (ptr+1):(ptr+n*n);
    yearVec(idx) = int32(yr);
    iVec(idx)    = int16(I);
    jVec(idx)    = int16(J);
    aVec(idx)    = vecA;
    zVec(idx)    = vecZ;
    xjVec(idx)   = vecXj;

    ptr = ptr + n*n;
end

A_panel = table(yearVec, iVec, jVec, aVec, zVec, xjVec, ...
    'VariableNames', {'year','i','j','a','z_ij','x_j'});

outFile = fullfile(outDir,'A_panel.csv');
writetable(A_panel, outFile);

fprintf('Done. Wrote: %s\n', outFile);
fprintf('Also wrote: %s\n', fullfile(outDir,'sector_map.csv'));
