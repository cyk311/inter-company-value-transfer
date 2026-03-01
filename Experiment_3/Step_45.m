%% 第四第五部分
% Step 4-5 (MATLAB) with Perron-Frobenius closure:
% 1) lambda = l*(I-A)^(-1)  (unchanged)
% 2) Solve p and (1+r) from p = p*(A + w*l)*(1+r)
%    => p*(A + w*l) = mu * p, mu = 1/(1+r)
%    p is the left PF eigenvector of B = A + w*l
% 3) w from IO: resident consumption subtotal column Y9:Y26 / Y27
% 4) Normalize using total price = total value: p*x = lambda*x
% 5) Output panel: year,i,j,c_ij,a,a_hat,p_i,p_ij,lambda_i,ue_ij,r

clear; clc;

%% ===== USER SETTINGS =====
ioFile   = 'IO(1981-2018).xlsx';
predFile = 'Ahat_panel.csv';     % from Stata
outDir   = '.';

% Numeric block definition consistent with your "人大IO表" memory
rangeNum = 'D9:AI33';            % 25x32

yearsWanted = (1981:2018)';
n = 18;

% Row indices within numeric block (D9:AI33):
% 1-18 intermediate input rows
% 19 intermediate input total
% 20 labor compensation
% 25 total input
row_interTotal = 19;
row_laborComp  = 20;

% Column indices within numeric block:
% 1-18 intermediate use columns
% 32 total output column (AI)
col_totalOut   = 32;

% "居民消费小计" is in Excel column Y per your instruction.
% Within D..AI range, Y maps to index: 25-4+1 = 22
col_resConsSub = 22;
%% ======== 模型初始化 ========

assert(isfile(ioFile),   'Cannot find %s', ioFile);
assert(isfile(predFile), 'Cannot find %s', predFile);

pred = readtable(predFile);
% Expect columns: year,i,j,a,a_hat,c_ij,model (model optional)
needVars = {'year','i','j','a','a_hat','c_ij'};
for k = 1:numel(needVars)
    assert(any(strcmp(pred.Properties.VariableNames, needVars{k})), ...
        'Ahat_panel.csv must contain column: %s', needVars{k});
end

pred.year = double(pred.year);
pred.i    = double(pred.i);
pred.j    = double(pred.j);

years = intersect(yearsWanted, unique(pred.year));
assert(~isempty(years), 'No overlapping years between predictions and 1981-2018.');

allSheets = sheetnames(ioFile);

% Preallocate outputs (same row count as pred)
Nobs = height(pred);
p_i      = nan(Nobs,1);
p_ij     = nan(Nobs,1);
lambda_i = nan(Nobs,1);
ue_ij    = nan(Nobs,1);
r_vec    = nan(Nobs,1);

% Sort by year for efficient processing
pred.idx0 = (1:Nobs)';
pred = sortrows(pred, {'year','j','i'}); % j,i order matches reshape column-major
p_i      = p_i(pred.idx0);   % dummy placeholders; will reassign later
p_ij     = p_ij(pred.idx0);
lambda_i = lambda_i(pred.idx0);
ue_ij    = ue_ij(pred.idx0);
r_vec    = r_vec(pred.idx0);

% We will fill in the sorted order then unsort at end.
out_p_i      = nan(Nobs,1);
out_p_ij     = nan(Nobs,1);
out_lambda_i = nan(Nobs,1);
out_ue_ij    = nan(Nobs,1);
out_r        = nan(Nobs,1);

fprintf('Processing years: %d to %d (%d years)\n', years(1), years(end), numel(years));

%% ======== 主循环 ========

for tt = 1:numel(years)
    yr = years(tt);
    sh = num2str(yr);

    if ~any(strcmp(allSheets, sh))
        warning('Sheet %s not found, skipping year %d.', sh, yr);
        continue
    end

    % Read numeric block
    M = readmatrix(ioFile,'Sheet',sh,'Range',rangeNum);

    % A from intermediate transactions Z and gross output x
    Z = M(1:n, 1:n);
    x = M(1:n, col_totalOut);          % gross output (n x 1)
    if any(~isfinite(x)) || any(x<=0)
        warning('Year %d: invalid gross output x; skipping.', yr);
        continue
    end

    % 构造技术矩阵
    A = Z ./ (x');                     % A_ij = z_ij / x_j
    
    % l for lambda: keep your prior method (labor compensation / output) (UNCHANGED)
    labComp = M(row_laborComp, 1:n)';  % (n x 1)
    l_row = (labComp ./ x)';           % (1 x n)

    % 计算价值向量
    % lambda = l*(I-A)^(-1)
    IA = eye(n) - A;
    lambda_row = l_row / IA;           % (1 x n)
    lambda_col = lambda_row(:);        % (n x 1)

    % 工人消费向量的构建使用居民消费支出列的结构
    % w from resident consumption subtotal column:
    % w = (Y9:Y26)/(Y27) -> within block: rows 1:18 over row 19, column col_resConsSub
    w_num = M(1:n, col_resConsSub);
    w_den = M(row_interTotal, col_resConsSub);
    if ~isfinite(w_den) || w_den==0 || any(~isfinite(w_num))
        warning('Year %d: invalid resident-consumption subtotal column; using uniform w.', yr);
        w = ones(n,1) / n;
    else
        w = w_num ./ w_den;            % (n x 1)
        % basic sanity: enforce nonnegative and sum to 1 (numerical tidy)
        w(w<0) = 0;
        s = sum(w);
        if s<=0
            w = ones(n,1) / n;
        else
            w = w ./ s;
        end
    end
    

    % Build B = A + w*l  (outer product)
    % p=p(1+r)(A+fl)=p(1+r)B
    B = A + (w * l_row);               % (n x n)

    % Perron-Frobenius: left eigenvector p of B (equivalently right eigenvector of B')
    [V,D] = eig(B');                   % B' * v = mu * v
    evals = diag(D);
    % 左特征向量即为p

    % Choose eigenvalue with largest real part (PF should be real & dominant for positive B)
    [~, idxMax] = max(real(evals));
    mu = real(evals(idxMax));
    v  = V(:, idxMax);

    % Convert to real vector; handle small imaginary parts
    if norm(imag(v)) > 1e-8*max(1,norm(real(v)))
        warning('Year %d: PF eigenvector has non-negligible imaginary part; taking real().', yr);
    end
    p_col = real(v);         %只留实部作为价格向量
    % Fix sign: make it mostly positive
    if sum(p_col) < 0        %价格之和如果小于0
        p_col = -p_col;      %反转价格向量元素的符号
    end
    % If still has negative entries due to reducibility / numerical issues, floor tiny negatives
    p_col(p_col < 0 & abs(p_col) < 1e-10) = 0;       %对绝对值特别小的负价格做归零处理
    
    %检查特征值和特征向量有效性，特征值 mu=1/(1+r) 为正，p不能为0
    if mu<=0 || ~isfinite(mu) || all(p_col==0)
        warning('Year %d: invalid PF eigenvalue/vector; skipping.', yr);
        continue
    end

    % From p = p*B*(1+r): pB = mu p, mu = 1/(1+r) => r = 1/mu - 1
    r = (1/mu) - 1;

    % Normalize p scale by total price = total value:
    % total value = lambda*x, total price = p*x
    p_row_raw = p_col';                % 1 x n
    totalValue = lambda_row * x;       % scalar
    totalPrice = p_row_raw  * x;       % scalar

    if ~isfinite(totalPrice) || totalPrice==0
        warning('Year %d: totalPrice invalid; skipping.', yr);
        continue
    end
    kScale = totalValue / totalPrice;  % 缩放因子，使得总价值=总价格
    p_row = p_row_raw * kScale;        % normalized p (1 x n)
    p_col = p_row(:);                  % (n x 1)

    %% 计算p_ij和UE_ij
    % Pull this year's panel rows (pred is sorted by year,j,i)
    idxYear = find(pred.year==yr);
    if isempty(idxYear)
        continue
    end

    Ty = pred(idxYear, :);
    % Ensure complete n*n block; if not, we compute row-wise without reshape
    % 检查第一象限是不是方阵
    % 如果不是方阵，则：
    if height(Ty) ~= n*n
        % Row-wise (robust fallback)
        for rr = 1:height(Ty)
            ii = Ty.i(rr); jj = Ty.j(rr);
            a0 = Ty.a(rr); ah = Ty.a_hat(rr);

            out_p_i(idxYear(rr))      = p_col(ii);
            out_lambda_i(idxYear(rr)) = lambda_col(ii);
            out_r(idxYear(rr))        = r;

            if isfinite(a0) && a0~=0 && isfinite(ah)
                pij = p_col(ii) * (ah / a0);            % 计算p_ij：用 â/a 修正 p_i 得到 i←j 的“交换价格”
                out_p_ij(idxYear(rr)) = pij;            % 输出 p_ij
                out_ue_ij(idxYear(rr)) = pij / lambda_col(ii);     % 输出 UE_ij
            end
        end
    
    % 如果第一象限是方阵，则：
    else
        % Matrix way (fast): reshape in (i,j) order
        aMat    = reshape(Ty.a,     n, n);   % i x j
        ahatMat = reshape(Ty.a_hat, n, n);

        ratio = nan(n,n);
        mask = isfinite(aMat) & (aMat~=0) & isfinite(ahatMat);
        ratio(mask) = ahatMat(mask) ./ aMat(mask);

        % p_i repeated across columns j
        pMat = repmat(p_col, 1, n);          % n x n
        lamMat = repmat(lambda_col, 1, n);   % n x n
        
        % 构造双边交换价格矩阵和双边不平等交换矩阵
        pijMat = pMat .* ratio;
        ueMat  = pijMat ./ lamMat;
        
        % 输出
        out_p_i(idxYear)      = pMat(:);
        out_lambda_i(idxYear) = lamMat(:);
        out_p_ij(idxYear)     = pijMat(:);
        out_ue_ij(idxYear)    = ueMat(:);
        out_r(idxYear)        = r;
    end

    fprintf('Year %d done. r=%.6g, mu=%.6g\n', yr, r, mu);
end

% Assemble final panel output in requested schema
out = table();
out.year     = pred.year;
out.i        = pred.i;
out.j        = pred.j;
out.c_ij     = pred.c_ij;
out.a        = pred.a;
out.a_hat    = pred.a_hat;
out.p_i      = out_p_i;
out.p_ij     = out_p_ij;
out.lambda_i = out_lambda_i;
out.ue_ij    = out_ue_ij;
out.r        = out_r;

% Restore original ordering (by original row order of Ahat_panel.csv)
% pred.idx0 stored original indices; but we overwrote pred.idx0 earlier—rebuild properly:
% We can reconstruct by re-reading and merging on (year,i,j,c_ij,a,a_hat) is risky.
% Instead, we kept only sorted pred. We'll output sorted by year,i,j which is typically OK for panel.
out = sortrows(out, {'year','i','j'});

outFile = fullfile(outDir, 'value_transfer_panel_pf.csv');
writetable(out, outFile);
fprintf('Done. Wrote: %s\n', outFile);
