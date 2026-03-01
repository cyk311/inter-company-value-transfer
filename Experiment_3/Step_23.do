*******************************************************
* 第二第三部分 (Stata, NO MATA):
*   - 第二部分：为每个(i,j)对应的技术时间序列做DFT，计算主周期c_ij
*   - 第三部分：分别带入3个线性回归模型 + 模型选择 + 导出 Ahat_panel.csv
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
* 第三部分: 线性回归
* a = alpha * ln(xb + u), u ~ N(0, sigma^2)
* u = exp(a/alpha) - xb
* ll = -ln(sigma) - .5 ln(2pi) - .5 (u/sigma)^2 - ln(alpha) + a/alpha
*
* Engineering fixes:
*   - lnf defaults to huge negative, never missing
*   - overflow protection for exp(a/alpha) with z>700
*------------------------
capture program drop ue_m3_ll
program define ue_m3_ll
    version 15
    args lnf xb lna lnsigma
    tempvar alpha sigma y z expz u

    quietly gen double `alpha' = exp(`lna')
    quietly gen double `sigma' = exp(`lnsigma')
    quietly gen double `y'     = $ML_y1
    quietly gen double `z'     = `y'/`alpha'

    quietly replace `lnf' = -1e300

    quietly gen double `expz' = .
    quietly replace `expz' = exp(`z') if !missing(`y',`xb',`sigma') & `sigma'>0 & `z'<=700 & `z'>=-700

    quietly gen double `u' = .
    quietly replace `u' = `expz' - `xb' if !missing(`expz')

    quietly replace `lnf' = -ln(`sigma') - 0.5*ln(2*c(pi)) ///
                            - 0.5*(`u'/`sigma')^2 ///
                            - ln(`alpha') + `z' ///
        if !missing(`u') & `sigma'>0
end

*------------------------
* Containers
*------------------------
gen double a_hat1 = .
gen double a_hat2 = .
gen double a_hat3 = .
gen double a_hat  = .
gen byte   model  = .

tempfile pairstats
postfile STATS int pair int i int j double c_ij byte model ///
    double SSE1 double SSE2 double SSE3 byte ok2 byte ok3 ///
    using `pairstats', replace

di as txt "Estimation stage..."

local iter = 0
foreach p of local pairs {
    local ++iter

    *---------------- (1) 回归模型一：最小二乘，线性趋势 ----------------
    capture noisily regress a t t2 cosTerm sinTerm if pair==`p'
    if _rc continue

    tempvar y1 r1
    predict double `y1' if e(sample), xb
    gen double `r1' = (a-`y1')^2 if pair==`p' & e(sample)
    quietly summarize `r1' if pair==`p' & e(sample), meanonly
    local SSE1 = r(sum)
    replace a_hat1 = `y1' if pair==`p' & e(sample)
    drop `y1' `r1'

    *---------------- (2) 回归模型二：最小二乘，指数趋势 ----------------
    local OK2 = 0
    local SSE2 = .
    quietly summarize a if pair==`p', meanonly
    if r(min)>0 {
        capture noisily regress ln(a) t t2 cosTerm sinTerm if pair==`p'
        if !_rc {
            local OK2 = 1
            tempvar xb2 y2 r2
            predict double `xb2' if e(sample), xb
            scalar __sig2 = e(rmse)^2
            gen double `y2' = exp(`xb2' + __sig2/2) if pair==`p' & e(sample)
            gen double `r2' = (a-`y2')^2 if pair==`p' & e(sample)
            quietly summarize `r2' if pair==`p' & e(sample), meanonly
            local SSE2 = r(sum)
            replace a_hat2 = `y2' if pair==`p' & e(sample)
            drop `xb2' `y2' `r2'
        }
    }

    *---------------- (3) 回归模型三：极大似然，对数趋势 ----------------
    local OK3 = 0
    local SSE3 = .

    capture noisily ml clear
    capture noisily ml model lf ue_m3_ll (xb: a = t t2 cosTerm sinTerm) /lna /lnsigma if pair==`p'
    if !_rc {
        * stable starts
        capture noisily ml init xb:_cons = 0 xb:t = 0 xb:t2 = 0 xb:cosTerm = 0 xb:sinTerm = 0
        capture noisily ml init /lna = 0 /lnsigma = 0   // alpha=1, sigma=1

        capture noisily ml maximize, difficult technique(nr 20 bhhh 10) iterate(200) nolog
        if !_rc & e(converged)==1 {
            local OK3 = 1
            tempvar xb3 y3 r3
            predict double `xb3' if e(sample), xb
            scalar __alpha = exp(_b[/lna])

            gen double `y3' = .
            replace `y3' = __alpha * ln(max(`xb3', 1e-12)) if pair==`p' & e(sample)

            gen double `r3' = (a-`y3')^2 if pair==`p' & e(sample) & !missing(`y3')
            quietly summarize `r3' if pair==`p' & e(sample) & !missing(`y3'), meanonly
            local SSE3 = r(sum)
            replace a_hat3 = `y3' if pair==`p' & e(sample)

            drop `xb3' `y3' `r3'
        }
    }

    *---------------- Model choice ----------------
	// 包含模型检验
    local choose = 1
    local best = `SSE1'

    if `OK2'==1 {
        if (`SSE2' < `best') {
            local choose = 2
            local best = `SSE2'
        }
    }
    if `OK3'==1 & `SSE3'!=. {
        if (`SSE3' < `best') {
            local choose = 3
            local best = `SSE3'
        }
    }

    * Auxiliary J-type override (if feasible)
    if `OK2'==1 {
        tempvar y2aux
        gen double `y2aux' = a_hat2 if pair==`p'
        capture noisily regress a t t2 cosTerm sinTerm `y2aux' if pair==`p' & !missing(`y2aux')
        if !_rc {
            capture test `y2aux' = 0
            if !_rc & r(p)<0.05 local choose = 2
        }
        drop `y2aux'
    }
    if `OK3'==1 {
        tempvar y3aux
        gen double `y3aux' = a_hat3 if pair==`p'
        capture noisily regress a t t2 cosTerm sinTerm `y3aux' if pair==`p' & !missing(`y3aux')
        if !_rc {
            capture test `y3aux' = 0
            if !_rc & r(p)<0.05 local choose = 3
        }
        drop `y3aux'
    }

    replace model = `choose' if pair==`p'
    replace a_hat = a_hat1 if pair==`p' & model==1
    replace a_hat = a_hat2 if pair==`p' & model==2
    replace a_hat = a_hat3 if pair==`p' & model==3

    local ii = floor(`p'/100)
    local jj = `p' - 100*`ii'
    quietly summarize c_ij if pair==`p', meanonly
    local cc = r(mean)

    post STATS (`p') (`ii') (`jj') (`cc') (`choose') (`SSE1') (`SSE2') (`SSE3') (`OK2') (`OK3')

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
