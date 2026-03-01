%% run_T_Rho_PF.m
%  目标：
%   - 年度 T 矩阵：T_{i<-j}^t = Z_{ij}^t * (1 - 1/rho_i^t),  rho_i^t = p_i^t / lambda_i^t
%   - 输出：
%       1) T_matrices_1981_2018.xlsx : 每年一张sheet(18×18)
%       2) T_panel_1981_2018.xlsx    : 面板长表 Year,i,j,T
%       3) rho_matrix_1981_2018.xlsx : rho矩阵 (Years×18)
%
%  关键计算口径：
%   (A) 进出口口径A：出口AE9:AE26、进口AG9:AG26 外生读取，但不用于A构造（仅保留读取接口）
%   (B) 价值型IO修正（冯志轩/ Marelli）且忽略折旧D：
%       A_{ij} = Z_{ij} / x_j
%       A*_{ij}= A_{ij} * (m_j/m_i)
%       l*_{i} = 劳动者报酬_i / m_i
%       lambda = l* (I - A*)^{-1}
%   (C) 生产价格与平均利润率：Perron–Frobenius
%       p A* = mu p  => 1+r = 1/mu
%       再用总价格=总价值归一化：p m = lambda m

clear; clc;

%% ========== 用户配置 ==========
inputFile = 'IO(1981-2018).xlsx';
years     = 1981:2018;
N         = 18;

% 数值区（25×32）
rangeNumeric = 'D9:AI33';

% 进出口（外生）
rangeExport  = 'AE9:AE26';
rangeImport  = 'AG9:AG26';

% 输出文件
outTmat  = 'T_matrices_1981_2018.xlsx';
outTpan  = 'T_panel_1981_2018.xlsx';
outRhoMx = 'rho_matrix_1981_2018.xlsx';

% 数值稳定
rcondThreshold = 1e-12;
epsDenom       = 1e-12;

%% ========== 清理旧输出 ==========
safeDelete(outTmat);
safeDelete(outTpan);
safeDelete(outRhoMx);

%% ========== 面板T预分配 ==========
Tcount = numel(years) * N * N;
YearVec = zeros(Tcount,1);
iVec    = zeros(Tcount,1);
jVec    = zeros(Tcount,1);
Tvec    = zeros(Tcount,1);
lambdaVec = zeros(Tcount,1);
pVec    = zeros(Tcount,1);

[II, JJ] = ndgrid(1:N, 1:N);
II = II(:); JJ = JJ(:);

ptr = 0;

%% ========== rho矩阵 ==========
rhoMat = zeros(numel(years), N);

%% ========== 主循环 ==========
for yy = 1:numel(years)
    t = years(yy);
    sheetName = num2str(t);

    % ---- 读取数值区 ----
    M = readmatrix(inputFile, 'Sheet', sheetName, 'Range', rangeNumeric);
    M = sanitizeNumeric(M);

    % ---- 读取外生进出口（口径A：不参与计算，仅保留接口） ----
    exportVec = tryReadVector(inputFile, sheetName, rangeExport, N); %#ok<NASGU>
    importVec = tryReadVector(inputFile, sheetName, rangeImport, N); %#ok<NASGU>

    % ---- 行列定位（人大IO表结构） ----
    ROW_INT_SUM = 19;
    ROW_COMP    = 20;
    ROW_TOTIN   = 25;
    COL_TOTOUT  = 32;   % 行总产出

    % ---- 提取基础量 ----
    Z = M(1:N, 1:N);                 % 18×18 中间投入流量（价值型）
    x_col = M(ROW_TOTIN, 1:N);       % 列总投入（作为A分母）
    x_col = max(x_col, epsDenom);

    % m：部门产品“总市场价格量”（用行总产出代理）
    m = M(1:N, COL_TOTOUT);
    m = max(m(:), epsDenom);

    comp   = M(ROW_COMP, 1:N);       % 劳动者报酬（按部门列）
    intSum = M(ROW_INT_SUM, 1:N);    % 中间投入合计（按部门列） %#ok<NASGU>

    %% ===== 1) 构造 A（价值型）并做 Marelli/Feng 修正得到 A*、l* =====
    A = Z ./ x_col;                  % A_{ij} = Z_{ij}/x_j

    ratioMat = (m(:)') ./ m(:);      % (m_j/m_i)
    Astar = A .* ratioMat;           % A*_{ij} = A_{ij}*(m_j/m_i)

    lstar = comp(:)' ./ m(:)';       % l*_i = 劳动者报酬_i / m_i（忽略折旧）

    %% ===== 2) 价值：lambda = l* (I - A*)^{-1} =====
    I = eye(N);
    B = I - Astar;
    if rcond(B) < rcondThreshold
        lambda = lstar * pinv(B);
        warning('Year %d: (I-A*) 病态，使用 pinv。rcond=%.3e', t, rcond(B));
    else
        lambda = lstar / B;
    end
    lambda = max(lambda, epsDenom);  % 防止后续除零

    %% ===== 3) 生产价格 & 平均利润率：Perron–Frobenius =====
    % 求 A* 的左Perron向量 p_raw（行向量）与主特征值 mu
    [mu, p_raw] = leftPerron(Astar);

    % 平均利润率：1+r = 1/mu
    r = (1/mu) - 1; %#ok<NASGU>  % 如果你想输出r，可在此处保存

    % 用“总价格=总价值”归一化：p m = lambda m
    valueTotal = lambda * m;         % 标量
    priceTotal = p_raw * m;          % 标量
    if abs(priceTotal) < epsDenom
        warning('Year %d: p_raw*m 接近0，改用单位缩放（p_raw不变）', t);
        scale = 1;
    else
        scale = valueTotal / priceTotal;
    end
    p = scale * p_raw;               % 归一化后的生产价格向量（1×N）
    p = max(p, epsDenom);

    %% ===== 4) rho 与 T =====
    rho = p ./ lambda;               % 1×N
    rhoMat(yy, :) = rho;

    % T_{i<-j} = Z_{ij} * (1 - 1/rho_i)（按行i缩放）
    factor = exp(p) ./ exp(lambda);
    T = Z .* factor(:);

    %% ===== 5) 输出：年度矩阵版T =====
    writematrix(T, outTmat, 'Sheet', sheetName, 'Range', 'A1');
    outLambdaMx = 'lambda_matrix_1981_2018.xlsx';
    outPMx      = 'p_matrix_1981_2018.xlsx';


    %% ===== 6) 输出：面板长表T =====
    idx = ptr + (1:(N*N));
    YearVec(idx) = t;
    iVec(idx)    = II;
    jVec(idx)    = JJ;
    Tvec(idx)    = T(:);
    lambdaVec(idx)  = lambda(:);
    pVec(idx)    = p(:);
    ptr = ptr + N*N;
end

%% ========== 输出面板T ==========
Tpanel = table(YearVec, iVec, jVec, lambdaVec, pVec, 'VariableNames', {'Year','i','j','lambda','p'});
writetable(Tpanel, outTpan, 'Sheet', 'panel', 'Range', 'A1');

%% ========== 输出rho矩阵（Years×18） ==========
% 形式：第一列Year，后18列rho_1...rho_18
rhoOut = [(years(:)), rhoMat];
header = [{'Year'}, arrayfun(@(k) sprintf('rho_%d', k), 1:18, 'UniformOutput', false)];
writecell(cellstr(header), outRhoMx, 'Sheet', 'rho', 'Range', 'A1');
writematrix(rhoOut, outRhoMx, 'Sheet', 'rho', 'Range', 'A2');

fprintf('\n完成。\nT矩阵版：%s\nT面板版：%s\nrho矩阵：%s\n', outTmat, outTpan, outRhoMx);

%% ====================== 本地函数 ======================
function safeDelete(fname)
    if exist(fname, 'file') == 2
        try
            delete(fname);
        catch
            warning('无法删除旧文件：%s（可能被Excel占用）。请关闭Excel后重试。', fname);
        end
    end
end

function X = sanitizeNumeric(X)
    if isempty(X)
        error('读取到空矩阵，请检查Range是否正确。');
    end
    X(isnan(X)) = 0;
end

function v = tryReadVector(file, sheet, range, N)
    try
        v = readmatrix(file, 'Sheet', sheet, 'Range', range);
        v = v(:);
        if numel(v) ~= N
            warning('读取向量维度异常：Sheet=%s Range=%s（期望%d，得到%d）', sheet, range, N, numel(v));
            v = nan(N,1);
        end
    catch
        warning('未能读取向量：Sheet=%s Range=%s（将置为NaN）', sheet, range);
        v = nan(N,1);
    end
end

function [mu, p_row] = leftPerron(A)
    % 返回非负矩阵A的主特征值 mu 与左Perron向量 p_row（行向量）
    % 通过对 A' 求主特征对：
    %   A' v = mu v  =>  v' A = mu v'
    try
        opts = struct();
        opts.tol = 1e-10;
        opts.maxit = 5000;
        [v, d] = eigs(A', 1, 'largestreal', opts);
        mu = real(d(1,1));
        v  = real(v);
    catch
        % eigs失败时回退到全特征分解
        [V, D] = eig(A');
        d = real(diag(D));
        [mu, k] = max(d);
        v = real(V(:,k));
    end

    % 规范化符号与非负性（PF向量可取同向非负）
    if sum(v) < 0
        v = -v;
    end
    v = abs(v);  % 防止数值误差导致的极小负数

    % 行向量形式
    p_row = (v(:))';
    % 若 mu 非正，强制报错（PF理论下非负矩阵主特征值应>=0）
    if ~(mu > 0)
        error('Perron特征值 mu 非正（mu=%.4g），请检查A*是否存在全零行/列或异常数据。', mu);
    end
end
