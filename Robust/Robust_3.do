*******************************************************
* Robust_3 (Stata, NO MATA)
* 任务：
* 1) 将时间分为 1981-1999 与 2000-2018 两段，分别按 DFT 提取周期
* 2) 使用同一回归方程：a = b0 + b1*t + b2*t2 + b3*cos + b4*sin
* 3) 每段都先在子样本区间估计，再对完整 1981-2018 生成预测值
* 4) 1981-1999 训练得到的完整期预测记作 ahat_1st_half
*    2000-2018 训练得到的完整期预测记作 ahat_2nd_half
* 5) 输出字段：year i j a ahat ahat_1st_half ahat_2nd_half
*******************************************************
version 15
clear all
set more off
set varabbrev off
set matsize 11000

***** USER SETTINGS *****
global INDIR  "E:\Bundle\2025学推\Robustness"
global OUTDIR "E:\Bundle\2025学推\Robustness"
local baseyear = 1981
*************************

*------------------------
* Load baseline panel from Step_23 output
*------------------------
import delimited using "${INDIR}/Ahat_panel.csv", clear varnames(1)
drop if missing(year,i,j)

gen int pair = i*100 + j
sort pair year

gen int t = year - `baseyear' + 1
gen double t2 = t^2

* Baseline ahat from Step_23 is assumed to be a_hat
capture confirm variable a_hat
if _rc {
    di as err "Variable a_hat not found in Ahat_panel.csv."
    exit 111
}

gen double ahat_1st_half = .
gen double ahat_2nd_half = .

* Segment-specific period containers
gen double c_1 = .
gen double c_2 = .

*======================================================
* Segment 1: Train on 1981-1999
*======================================================
local ystart = 1981
local yend   = 1999

levelsof pair if inrange(year,`ystart',`yend'), local(pairs_1)
local npairs_1 : word count `pairs_1'
di as txt "Segment 1 DFT: total pairs = `npairs_1'"

local iter = 0
foreach p of local pairs_1 {
    local ++iter

    preserve
        keep if pair==`p' & inrange(year,`ystart',`yend')
        sort year

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

        local K = floor(`N'/2)
        local bestmag = -1
        local bestk   = 1

        forvalues k = 1/`K' {
            gen double ang = 2*c(pi)*`k'*tt/`N'
            gen double cs  = cos(ang)
            gen double sn  = sin(ang)
            gen double xcs = x*cs
            gen double xsn = x*sn

            egen double Re = total(xcs)
            egen double Im = total(-xsn)

            quietly summarize Re in 1, meanonly
            local ReS = r(mean)
            quietly summarize Im in 1, meanonly
            local ImS = r(mean)

            local mag = sqrt((`ReS')^2 + (`ImS')^2)
            if (`mag' > `bestmag') {
                local bestmag = `mag'
                local bestk   = `k'
            }

            drop ang cs sn xcs xsn Re Im
        }

        local c = `N'/`bestk'
    restore

    replace c_1 = `c' if pair==`p'

    if mod(`iter',50)==0 di as txt "Segment 1 DFT progress: `iter' / `npairs_1'"
}

replace c_1 = 2 if missing(c_1) | c_1 < 2
gen double cos_1 = cos(2*c(pi)*t/c_1)
gen double sin_1 = sin(2*c(pi)*t/c_1)

local iter = 0
foreach p of local pairs_1 {
    local ++iter

    capture noisily regress a cos_1 sin_1 if pair==`p' & inrange(year,`ystart',`yend')
    if _rc continue

    tempvar yhat1
    predict double `yhat1' if pair==`p', xb
    replace ahat_1st_half = `yhat1' if pair==`p' & inrange(year,1981,2018)
    drop `yhat1'

    if mod(`iter',50)==0 di as txt "Segment 1 Est progress: `iter' / `npairs_1'"
}

*======================================================
* Segment 2: Train on 2000-2018
*======================================================
local ystart = 2000
local yend   = 2018

levelsof pair if inrange(year,`ystart',`yend'), local(pairs_2)
local npairs_2 : word count `pairs_2'
di as txt "Segment 2 DFT: total pairs = `npairs_2'"

local iter = 0
foreach p of local pairs_2 {
    local ++iter

    preserve
        keep if pair==`p' & inrange(year,`ystart',`yend')
        sort year

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

        local K = floor(`N'/2)
        local bestmag = -1
        local bestk   = 1

        forvalues k = 1/`K' {
            gen double ang = 2*c(pi)*`k'*tt/`N'
            gen double cs  = cos(ang)
            gen double sn  = sin(ang)
            gen double xcs = x*cs
            gen double xsn = x*sn

            egen double Re = total(xcs)
            egen double Im = total(-xsn)

            quietly summarize Re in 1, meanonly
            local ReS = r(mean)
            quietly summarize Im in 1, meanonly
            local ImS = r(mean)

            local mag = sqrt((`ReS')^2 + (`ImS')^2)
            if (`mag' > `bestmag') {
                local bestmag = `mag'
                local bestk   = `k'
            }

            drop ang cs sn xcs xsn Re Im
        }

        local c = `N'/`bestk'
    restore

    replace c_2 = `c' if pair==`p'

    if mod(`iter',50)==0 di as txt "Segment 2 DFT progress: `iter' / `npairs_2'"
}

replace c_2 = 2 if missing(c_2) | c_2 < 2
gen double cos_2 = cos(2*c(pi)*t/c_2)
gen double sin_2 = sin(2*c(pi)*t/c_2)

local iter = 0
foreach p of local pairs_2 {
    local ++iter

    capture noisily regress a cos_2 sin_2 if pair==`p' & inrange(year,`ystart',`yend')
    if _rc continue

    tempvar yhat2
    predict double `yhat2' if pair==`p', xb
    replace ahat_2nd_half = `yhat2' if pair==`p' & inrange(year,1981,2018)
    drop `yhat2'

    if mod(`iter',50)==0 di as txt "Segment 2 Est progress: `iter' / `npairs_2'"
}

*------------------------
* Hybrid
*------------------------

gen     ahat_hybrid = .
replace ahat_hybrid = ahat_1st_half if year >= 2000
replace ahat_hybrid = ahat_2nd_half if year <  2000

ttest 	 a_hat == ahat_hybrid
signrank a_hat =  ahat_hybrid

*------------------------
* 统计检验：ahat_hybrid vs a_hat
*------------------------
gen double diff_hat = ahat_hybrid - a_hat if !missing(ahat_hybrid,a_hat)

quietly count if !missing(diff_hat)
local Ndiff = r(N)

local tstat = .
local tpval = .
local mdiff = .
local sddiff = .

capture noisily ttest diff_hat == 0
if !_rc {
    local tstat  = r(t)
    local tpval  = r(p)
    local mdiff  = r(mu_1)
    local sddiff = r(sd_1)
}

local sr_z = .
local sr_p = .
capture noisily signrank diff_hat = 0
if !_rc {
    local sr_z = r(z)
    local sr_p = r(p)
}

*------------------------
* Export
*------------------------
	
preserve
    keep year i j a a_hat ahat_1st_half ahat_2nd_half ahat_hybrid
    rename a_hat ahat
    sort year i j

    gen str20 a_s          = string(a, "%20.15f")
    gen str20 ahat_s       = string(ahat, "%20.15f")
    gen str20 ahat1_s      = string(ahat_1st_half, "%20.15f")
    gen str20 ahat2_s      = string(ahat_2nd_half, "%20.15f")
	gen str20 ahath_s	   = string(ahat_hybrid, "%20.15f")
	
    keep year i j a_s ahat_s ahat1_s ahat2_s ahath_s
    rename a_s a
    rename ahat_s ahat
    rename ahat1_s ahat_1st_half
    rename ahat2_s ahat_2nd_half
	rename ahath_s ahat_hybrid

    export excel using "${OUTDIR}/Robust_3_compare.xlsx", ///
        sheet("data") firstrow(variables) replace
restore

putexcel set "${OUTDIR}/Robust_3_compare.xlsx", sheet("tests") modify
putexcel A1 = "test" B1 = "stat" C1 = "p_value" D1 = "N" E1 = "mean_diff" F1 = "sd_diff"
putexcel A2 = "paired_ttest" B2 = (`tstat') C2 = (`tpval') D2 = (`Ndiff') E2 = (`mdiff') F2 = (`sddiff')
putexcel A3 = "signrank" B3 = (`sr_z') C3 = (`sr_p') D3 = (`Ndiff')

di as txt "Done. Output written to:"
di as txt "  ${OUTDIR}/Robust_3_compare.xlsx"
*******************************************************
