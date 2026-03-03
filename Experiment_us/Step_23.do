*******************************************************
* 美国稳健性检验 Step_23 (优化版)
* 参考 Experiment_3/Step_23.do 的方法与回归模型：
*   1) 对每个(i,j)序列先去趋势后做DFT，得到主周期 c_ij
*   2) 回归模型：a_ij = alpha + beta1*t + beta2*t^2 + beta3*cos(2*pi*t/c_ij) + beta4*sin(2*pi*t/c_ij) + mu_t
*   3) 输出 Ahat_panel.csv 与 pair_diagnostics.csv
*
* 数据结构：year; i; j; a_ij
* 性能优化：
*   - 不使用 preserve/restore + 双重Stata循环
*   - 在 Mata 中一次性按 pair 分块计算 DFT 与回归拟合值
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
import delimited using "${INDIR}/A_panel.csv", clear varnames(1)

* 兼容 year; i; j; a_ij 结构
destring year i j a_ij, replace force
drop if missing(year, i, j)
rename a_ij a

* 使用 group() 构造高效 pair id（避免 i*100+j 可能冲突）
egen long pair = group(i j)
sort pair year

* time index
gen int t  = year - `baseyear' + 1
gen double t2 = t^2

* 预留结果变量
gen double c_ij  = .
gen double a_hat = .
gen byte   model = .

*------------------------
* Core compute in Mata (optimized for large panel)
*------------------------
mata:
real scalar has_missing(real colvector x) {
    return(sum(x :== .) > 0)
}

void step23_us_fast()
{
    real colvector pair = st_data(., "pair")
    real colvector a    = st_data(., "a")
    real colvector t    = st_data(., "t")
    real colvector t2   = st_data(., "t2")

    real scalar Nall
    real scalar s, e, N, K, k, bestk, n_pairs
    real scalar c, mag, bestmag

    real matrix X, Z
    real colvector y, trend, x, yh, re, im, SSE
    real colvector c_out, yhat_out, model_out

    Nall = rows(pair)
    c_out    = J(Nall, 1, 2)
    yhat_out = J(Nall, 1, .)
    model_out= J(Nall, 1, .)

    s = 1
    n_pairs = 0
    while (s <= Nall) {
        e = s
        while (e < Nall & pair[e+1] == pair[s]) e++
        n_pairs = n_pairs + 1

        y = a[|s\e|]
        N = rows(y)

        if (N >= 10 & !has_missing(y)) {
            X = J(N,1,1), t[|s\e|], t2[|s\e|]
            trend = X * qrsolve(X, y)
            x = y - trend

            K = floor(N/2)
            bestk = 1
            bestmag = -1

            for (k = 1; k <= K; k++) {
                real colvector ang
                ang = 2*pi()*k*(0::(N-1))/N
                re = x :* cos(ang)
                im = x :* sin(ang)
                mag = sqrt((sum(re))^2 + (sum(im))^2)
                if (mag > bestmag) {
                    bestmag = mag
                    bestk = k
                }
            }

            c = N / bestk
            if (c < 2 | c == .) c = 2

            c_out[|s\e|] = J(N,1,c)

            Z = J(N,1,1), t[|s\e|], t2[|s\e|], cos(2*pi()*t[|s\e|]/c), sin(2*pi()*t[|s\e|]/c)
            yh = Z * qrsolve(Z, y)

            yhat_out[|s\e|] = yh
            model_out[|s\e|] = J(N,1,1)
        }
        else {
            c_out[|s\e|] = J(N,1,2)
        }

        if (mod(n_pairs, 500) == 0) {
            printf("Processed pairs: %9.0f\r", n_pairs)
        }

        s = e + 1
    }

    st_store(., "c_ij", c_out)
    st_store(., "a_hat", yhat_out)
    st_store(., "model", model_out)

    printf("\nTotal pairs processed: %9.0f\n", n_pairs)
}

step23_us_fast()
end

replace c_ij = 2 if missing(c_ij) | c_ij < 2

*------------------------
* Diagnostics per pair
*------------------------
gen double sqe = (a - a_hat)^2 if model == 1

preserve
    collapse (first) i j c_ij model (sum) SSE1 = sqe, by(pair)
    order pair i j c_ij model SSE1
    export delimited using "${OUTDIR}/pair_diagnostics.csv", replace
restore

*------------------------
* Export panel outputs
*------------------------
keep year i j a a_hat c_ij model
order year i j a a_hat c_ij model
export delimited using "${OUTDIR}/Ahat_panel.csv", replace

di as txt "Done. Outputs written to:"
di as txt "  ${OUTDIR}/Ahat_panel.csv"
di as txt "  ${OUTDIR}/pair_diagnostics.csv"
*******************************************************
