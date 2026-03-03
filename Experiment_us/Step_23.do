*******************************************************
* 美国稳健性检验 Step_23 (高性能优化版)
* 参考 Experiment_3/Step_23.do 的模型：
*   1) 每个(i,j)序列去趋势后做DFT，得到主周期 c_ij
*   2) 回归: a_ij = alpha + beta1*t + beta2*t^2 + beta3*cos(2*pi*t/c_ij) + beta4*sin(2*pi*t/c_ij) + mu_t
* 数据结构: year; i; j; a_ij
*
* 关键优化：
*   - 用 Mata + panelsetup 一次性按 pair 分块处理（避免 preserve/restore 大循环）
*   - 用 FFT 代替逐频率显式求和，显著降低 DFT 计算量
*******************************************************
version 15
clear all
set more off
set varabbrev off
set matsize 11000

***** USER SETTINGS *****
global INDIR  "E:\Bundle\2025学推\Experiment_us"
global OUTDIR "E:\Bundle\2025学推\Experiment_us"
local baseyear = 1981
*************************

*------------------------
* Load panel
*------------------------
import delimited using "${INDIR}/A_panel.csv", clear varnames(1) stringcols(_all)
destring year i j a_ij, replace force
drop if missing(year,i,j)
rename a_ij a

* 稳健 pair id + 排序
egen long pair = group(i j)
sort pair year

* time index
gen int t  = year - `baseyear' + 1
gen double t2 = t^2

* 结果变量
gen double c_ij  = .
gen double a_hat = .
gen byte   model = .
gen double SSE1  = .

mata:
real scalar has_missing(real colvector x)
{
    return(sum(x :== .) > 0)
}

void step23_us_fft_fast()
{
    real colvector pair, a, t, t2
    real matrix P
    real scalar g, G, s, e, N, K, idx
    real scalar c, bestk, bestmag

    real matrix X, Z
    real colvector y, trend, x, yh, resid2, mag
    real colvector c_out, yhat_out, model_out, sse_out

    complex colvector F

    pair = st_data(., "pair")
    a    = st_data(., "a")
    t    = st_data(., "t")
    t2   = st_data(., "t2")

    P = panelsetup(pair, 1)
    G = rows(P)

    c_out     = J(rows(pair), 1, 2)
    yhat_out  = J(rows(pair), 1, .)
    model_out = J(rows(pair), 1, .)
    sse_out   = J(rows(pair), 1, .)


    for (g = 1; g <= G; g++) {
        s = P[g,1]
        e = P[g,2]
        N = e - s + 1

        y = a[|s\e|]

        if (N >= 10 & !has_missing(y)) {
            X = J(N,1,1), t[|s\e|], t2[|s\e|]
            trend = X * qrsolve(X, y)
            x = y - trend

            * FFT: F[1] 为零频; 频率 k 对应 F[k+1]
            F = fft(complex(x, J(N,1,0)))
            K = floor(N/2)

            if (K >= 1) {
                mag = abs(F[|2\(K+1)|])
                bestk = 1
                bestmag = mag[1]
                for (idx = 2; idx <= rows(mag); idx++) {
                    if (mag[idx] > bestmag) {
                        bestmag = mag[idx]
                        bestk = idx
                    }
                }
                if (rows(mag)==1) bestk = 1
            }
            else {
                bestk = 1
            }

            c = N / bestk
            if (c < 2 | c == .) c = 2
            c_out[|s\e|] = J(N,1,c)

            Z = J(N,1,1), t[|s\e|], t2[|s\e|], cos(2*pi()*t[|s\e|]/c), sin(2*pi()*t[|s\e|]/c)
            yh = Z * qrsolve(Z, y)
            yhat_out[|s\e|]  = yh
            model_out[|s\e|] = J(N,1,1)

            resid2 = (y - yh):^2
            sse_out[|s\e|] = J(N,1,sum(resid2))
        }
        else {
            c_out[|s\e|] = J(N,1,2)
        }

        if (mod(g, 1000)==0) {
            printf("Processed pairs: %9.0f\n", g)
        }
    }

    st_store(., "c_ij",  c_out)
    st_store(., "a_hat", yhat_out)
    st_store(., "model", model_out)
    st_store(., "SSE1",  sse_out)

    printf("Total pairs processed: %9.0f\n", G)
}

step23_us_fft_fast()
end

replace c_ij = 2 if missing(c_ij) | c_ij<2

*------------------------
* Export
*------------------------
preserve
    keep if model==1
    collapse (first) i j c_ij model SSE1, by(pair)
    order pair i j c_ij model SSE1
    export delimited using "${OUTDIR}/pair_diagnostics.csv", replace
restore

keep year i j a a_hat c_ij model
order year i j a a_hat c_ij model
export delimited using "${OUTDIR}/Ahat_panel.csv", replace

di as txt "Done. Outputs written to:"
di as txt "  ${OUTDIR}/Ahat_panel.csv"
di as txt "  ${OUTDIR}/pair_diagnostics.csv"
*******************************************************
