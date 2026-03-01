*******************************************************
* 第二第三部分 (Stata, NO MATA):
*   - 第二部分：为每个(i,j)对应的技术时间序列做DFT，计算主周期c_ij
*   - 第三部分：带入指定线性回归模型 + 导出 Ahat_panel.csv
*******************************************************
version 15
clear all
set more off
set varabbrev off
set matsize 11000

***** USER SETTINGS *****
global INDIR  "E:\Bundle\2025学推\Experiment_3"      // 存放 A_panel.csv 文件的文件夹位置(from MATLAB Step1)
global OUTDIR "E:\Bundle\2025学推\Experiment_3"      // 输出位置
local baseyear = 1981
*************************

*------------------------
* Load panel
*------------------------
import delimited using "${INDIR}/A_panel.csv", clear varnames(1) stringcols(_all)
destring year i j a z_ij x_j, replace force
drop if missing(year,i,j)

gen int pair = i*100 + j //构建(i,j)数组
sort pair year

* time index t = 1..T
gen int t  = year - `baseyear' + 1
gen double t2 = t^2

*------------------------
* 第二部分: 为每个(i,j)对应的技术时间序列做DFT，计算主周期c_ij (NO FFT, NO MATA)
* 对距平值做DFT： x_t = a_t - mean(a_t)
* For k=1..floor(N/2): magnitude = sqrt(Re^2 + Im^2)
*   实部：Re = sum_t x_t cos(2*pi*k*t/N), 虚部：Im = -sum_t x_t sin(2*pi*k*t/N)
* Choose k* with max magnitude; c = N/k*
*------------------------
gen double c_ij = .

levelsof pair, local(pairs)
local npairs : word count `pairs'
di as txt "DFT: total pairs = `npairs'"

local iter = 0    //迭代次数计数
foreach p of local pairs {  //对每个(i,j)数组循环
    local ++iter

    preserve
        keep if pair==`p'
        sort year

        * require no missing and enough obs
        count
        local N = r(N)
        quietly count if missing(a)
        local Nm = r(N)

        if (`N' < 10 | `Nm' > 0) {
            restore
            continue
        }

        quietly summarize a, meanonly
        gen double x = a - r(mean)
        gen int tt = _n - 1

        local K = floor(`N'/2) //需要遍历的周期的上下限
        local bestmag = -1
        local bestk   = 1

        * 循环遍历各个频率，找到强度最大的频率
        forvalues k = 1/`K' {
            gen double ang = 2*c(pi)*`k'*tt/`N'
            gen double cs  = cos(ang)
            gen double sn  = sin(ang)

            gen double xcs = x*cs
            gen double xsn = x*sn

            egen double Re = total(xcs)
            egen double Im = total(-xsn)

            * take the (constant) total from the first observation
            quietly summarize Re in 1, meanonly
            local ReS = r(mean)
            quietly summarize Im in 1, meanonly
            local ImS = r(mean)

            local mag = sqrt((`ReS')^2 + (`ImS')^2)
            if (`mag' > `bestmag') {
                local bestmag = `mag'
                local bestk   = `k'
            }
			* 记录最大的强度和对应的频率
			
            drop ang cs sn xcs xsn Re Im
			* 重置循环内的临时变量
        }

        local c = `N'/`bestk'
    restore

    replace c_ij = `c' if pair==`p'

    if mod(`iter',50)==0 di as txt "DFT progress: `iter' / `npairs'"
}

replace c_ij = 2 if missing(c_ij) | c_ij<2

* trig terms for Step 3
gen double cosTerm = cos(2*c(pi)*t/c_ij)
gen double sinTerm = sin(2*c(pi)*t/c_ij)

*------------------------
* 第三部分: 线性回归（仅保留模型一）
* a_{ij} = alpha + beta_1*t + beta_2*t^2 + beta_3*sin(2*pi*t/c_ij) + beta_4*cos(2*pi*t/c_ij) + mu_t
*------------------------
gen double a_hat  = .
gen byte   model  = .

tempfile pairstats
postfile STATS int pair int i int j double c_ij byte model ///
    double SSE1 ///
    using `pairstats', replace

di as txt "Estimation stage..."

local iter = 0
foreach p of local pairs {
    local ++iter

    *---------------- 回归模型：最小二乘，线性趋势 ----------------
    capture noisily regress a t t2 cosTerm sinTerm if pair==`p'
    if _rc continue

    tempvar y1 r1
    predict double `y1' if e(sample), xb
    gen double `r1' = (a-`y1')^2 if pair==`p' & e(sample)
    quietly summarize `r1' if pair==`p' & e(sample), meanonly
    local SSE1 = r(sum)
    replace a_hat = `y1' if pair==`p' & e(sample)
    drop `y1' `r1'

    replace model = 1 if pair==`p' & e(sample)

    local ii = floor(`p'/100)
    local jj = `p' - 100*`ii'
    quietly summarize c_ij if pair==`p', meanonly
    local cc = r(mean)

    post STATS (`p') (`ii') (`jj') (`cc') (1) (`SSE1')

    if mod(`iter',50)==0 di as txt "Est progress: `iter' / `npairs'"
}

postclose STATS

*------------------------
* Export
*------------------------
keep year i j a a_hat c_ij model
order year i j a a_hat c_ij model
export delimited using "${OUTDIR}/Ahat_panel.csv", replace

use `pairstats', clear
export delimited using "${OUTDIR}/pair_diagnostics.csv", replace

di as txt "Done. Outputs written to:"
di as txt "  ${OUTDIR}/Ahat_panel.csv"
di as txt "  ${OUTDIR}/pair_diagnostics.csv"
*******************************************************
