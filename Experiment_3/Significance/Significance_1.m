%% test_ue_direction_and_matrix.m
% Input : value_transfer_long.csv (columns include year,i,j,ue_ij)
% Output:
%   1) classMat:  18x18,  +1=inflow, -1=outflow, 0=none, NaN=insufficient
%   2) pMat:      18x18,  p-value (one-sided if directional, two-sided otherwise)
%   3) labelMat:  18x18 string, "+p", "-p", or "p", or "NA"
%   4) write to Excel: ue_significance_matrix.xlsx

clear; clc;

%% ===== User settings =====
inFile  = 'value_transfer_long.csv';
n       = 18;
alpha   = 0.05;   % significance level
minObs  = 10;     % minimum #years to run tests robustly
baseline = 1;     % UE baseline (no transfer direction)
outXlsx = 'ue_significance_matrix.xlsx';
% ==========================

assert(isfile(inFile), 'Cannot find %s', inFile);

T = readtable(inFile);

needVars = {'year','i','j','ue_ij'};
for k = 1:numel(needVars)
    assert(any(strcmp(T.Properties.VariableNames, needVars{k})), ...
        'Missing column %s in %s', needVars{k}, inFile);
end

% Ensure numeric
T.i = double(T.i);
T.j = double(T.j);
T.ue_ij = double(T.ue_ij);

% Containers
classMat = nan(n,n);          % +1, -1, 0, NaN
pMat     = nan(n,n);          % p-value
labelMat = strings(n,n);      % "+0.01", "-0.05", "0.23", "NA"
meanMat  = nan(n,n);          % mean(ue_ij)
nObsMat  = zeros(n,n);        % number of observations used

% Loop each (i,j)
for ii = 1:n
    for jj = 1:n
        idx = (T.i==ii) & (T.j==jj);
        x = T.ue_ij(idx);

        % Clean
        x = x(isfinite(x));
        nObsMat(ii,jj) = numel(x);

        if numel(x) < minObs
            classMat(ii,jj) = NaN;
            pMat(ii,jj) = NaN;
            labelMat(ii,jj) = "NA";
            continue
        end

        meanMat(ii,jj) = mean(x);

        d = x - baseline;  % test mean(d)=0

        % Handle degenerate variance (all constant)
        if all(abs(d - d(1)) < 1e-12)
            if d(1) > 0
                classMat(ii,jj) = +1;
                pMat(ii,jj) = 0;
                labelMat(ii,jj) = "+0";
            elseif d(1) < 0
                classMat(ii,jj) = -1;
                pMat(ii,jj) = 0;
                labelMat(ii,jj) = "-0";
            else
                classMat(ii,jj) = 0;
                pMat(ii,jj) = 1;
                labelMat(ii,jj) = "1";
            end
            continue
        end

        % One-sided tests for direction
        % H0: mean(d)=0
        % H1R: mean(d)>0  (net inflow)
        % H1L: mean(d)<0  (net outflow)
        [hR, pR] = ttest(d, 0, 'Alpha', alpha, 'Tail', 'right');
        [hL, pL] = ttest(d, 0, 'Alpha', alpha, 'Tail', 'left');

        if hR==1
            classMat(ii,jj) = +1;
            pMat(ii,jj) = pR;
            labelMat(ii,jj) = "+" + formatP(pR);
        elseif hL==1
            classMat(ii,jj) = -1;
            pMat(ii,jj) = pL;
            labelMat(ii,jj) = "-" + formatP(pL);
        else
            % No significant direction: report two-sided p-value
            p2 = min(1, 2*min(pR, pL));
            classMat(ii,jj) = 0;
            pMat(ii,jj) = p2;
            labelMat(ii,jj) = formatP(p2);
        end
    end
end

% Display quick summary
fprintf('Alpha = %.3f, minObs = %d\n', alpha, minObs);
fprintf('Inflow count  : %d\n', sum(classMat(:)==+1, 'omitnan'));
fprintf('Outflow count : %d\n', sum(classMat(:)==-1, 'omitnan'));
fprintf('No-dir count  : %d\n', sum(classMat(:)==0,  'omitnan'));
fprintf('NA count      : %d\n', sum(isnan(classMat(:))));

% ===== Write outputs (string-safe for older MATLAB) =====

% Build label sheet as cell array (n+1 by n+1)
hdr = cell(1, n+1);
hdr{1} = '';  % top-left blank
for k = 1:n
    hdr{1+k} = num2str(k);
end

outCell = cell(n+1, n+1);
outCell(1,:) = hdr;

for ii = 1:n
    outCell{1+ii,1} = num2str(ii);
    for jj = 1:n
        % labelMat is string array -> convert to char
        outCell{1+ii,1+jj} = char(labelMat(ii,jj));
    end
end

% Write label matrix
writecell(outCell, outXlsx, 'Sheet', 'label');

% -------- Numeric sheets: use writematrix + writecell for headers/row labels --------

% Helper headers/row labels (cell of char)
colHdr = cell(1,n); for k=1:n, colHdr{k}=num2str(k); end
rowHdr = cell(n,1); for k=1:n, rowHdr{k}=num2str(k); end

% p_value
writecell(colHdr, outXlsx, 'Sheet', 'p_value', 'Range', 'B1');
writecell(rowHdr, outXlsx, 'Sheet', 'p_value', 'Range', 'A2');
writematrix(pMat,   outXlsx, 'Sheet', 'p_value', 'Range', 'B2');

% class
writecell(colHdr, outXlsx, 'Sheet', 'class', 'Range', 'B1');
writecell(rowHdr, outXlsx, 'Sheet', 'class', 'Range', 'A2');
writematrix(classMat, outXlsx, 'Sheet', 'class', 'Range', 'B2');

% mean_ue
writecell(colHdr, outXlsx, 'Sheet', 'mean_ue', 'Range', 'B1');
writecell(rowHdr, outXlsx, 'Sheet', 'mean_ue', 'Range', 'A2');
writematrix(meanMat, outXlsx, 'Sheet', 'mean_ue', 'Range', 'B2');

% n_obs
writecell(colHdr, outXlsx, 'Sheet', 'n_obs', 'Range', 'B1');
writecell(rowHdr, outXlsx, 'Sheet', 'n_obs', 'Range', 'A2');
writematrix(nObsMat, outXlsx, 'Sheet', 'n_obs', 'Range', 'B2');

fprintf('Wrote: %s\n', outXlsx);


%% ===== helper: format p-values =====
function s = formatP(p)
    if isnan(p)
        s = "NA";
    elseif p==0
        s = "0";
    elseif p < 1e-4
        s = string(sprintf('%.2e', p));
    else
        s = string(sprintf('%.4f', p));
    end
end
