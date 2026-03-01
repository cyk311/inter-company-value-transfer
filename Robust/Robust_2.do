*******************************************************
* Robust_2 (Stata, NO MATA)
* 任务：
* 1) 将时间拆分为 1981-1999 与 2000-2018，两段分别按 Step_23.do 方法估计 a_hat_2
* 2) 比较 a_hat_2 与 Step_23.do 的 a_hat 是否存在统计差异
* 3) 输出 year, a, a_hat, a_hat_2 到同一张 xlsx
*******************************************************
version 15
clear all
set more off
set varabbrev off
set matsize 11000

***** USER SETTINGS *****
global INDIR  "E:\Bundle\2025学推\Experiment_3"
global OUTDIR "E:\Bundle\2025学推\Robust"
local baseyear = 1981
*************************

*------------------------
* Load panel
*------------------------
import delimited using "${INDIR}/A_panel.csv", clear varnames(1) stringcols(_all)
destring year i j a z_ij x_j, replace force
drop if missing(year,i,j)

gen int pair = i*100 + j
sort pair year

gen int t  = year - `baseyear' + 1
gen double t2 = t^2

sort year i j
* Merge Step_23 baseline a_hat
merge 1:1 year i j using "${INDIR}/Ahat_panel.csv", keepusing(a_hat) nogen

* Robust estimate container
gen double a_hat_2 = .
gen double c_ij_2  = .

*======================================================
* 分段估计：1981-1999
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

    replace c_ij_2 = `c' if pair==`p' & inrange(year,`ystart',`yend')

    if mod(`iter',50)==0 di as txt "Segment 1 DFT progress: `iter' / `npairs_1'"
}

replace c_ij_2 = 2 if inrange(year,`ystart',`yend') & (missing(c_ij_2) | c_ij_2<2)

gen double cosTerm2 = .
gen double sinTerm2 = .
replace cosTerm2 = cos(2*c(pi)*t/c_ij_2) if inrange(year,`ystart',`yend')
replace sinTerm2 = sin(2*c(pi)*t/c_ij_2) if inrange(year,`ystart',`yend')

local iter = 0
foreach p of local pairs_1 {
    local ++iter
    capture noisily regress a t t2 cosTerm2 sinTerm2 if pair==`p' & inrange(year,`ystart',`yend')
    if _rc continue

    tempvar y1
    predict double `y1' if e(sample), xb
    replace a_hat_2 = `y1' if pair==`p' & inrange(year,`ystart',`yend') & e(sample)
    drop `y1'

    if mod(`iter',50)==0 di as txt "Segment 1 Est progress: `iter' / `npairs_1'"
}

*======================================================
* 分段估计：2000-2018
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

    replace c_ij_2 = `c' if pair==`p' & inrange(year,`ystart',`yend')

    if mod(`iter',50)==0 di as txt "Segment 2 DFT progress: `iter' / `npairs_2'"
}

replace c_ij_2 = 2 if inrange(year,`ystart',`yend') & (missing(c_ij_2) | c_ij_2<2)

replace cosTerm2 = cos(2*c(pi)*t/c_ij_2) if inrange(year,`ystart',`yend')
replace sinTerm2 = sin(2*c(pi)*t/c_ij_2) if inrange(year,`ystart',`yend')

local iter = 0
foreach p of local pairs_2 {
    local ++iter
    capture noisily regress a t t2 cosTerm2 sinTerm2 if pair==`p' & inrange(year,`ystart',`yend')
    if _rc continue

    tempvar y2
    predict double `y2' if e(sample), xb
    replace a_hat_2 = `y2' if pair==`p' & inrange(year,`ystart',`yend') & e(sample)
    drop `y2'

    if mod(`iter',50)==0 di as txt "Segment 2 Est progress: `iter' / `npairs_2'"
}

*------------------------
* 统计检验：a_hat_2 vs a_hat
*------------------------
gen double diff_hat = a_hat_2 - a_hat if !missing(a_hat_2,a_hat)

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
* 导出 Excel
*------------------------
preserve
    keep year a a_hat a_hat_2
    sort year
    export excel using "${OUTDIR}/Robust_2_compare.xlsx", ///
        sheet("data") firstrow(variables) replace
restore

putexcel set "${OUTDIR}/Robust_2_compare.xlsx", sheet("tests") modify
putexcel A1 = "test" B1 = "stat" C1 = "p_value" D1 = "N" E1 = "mean_diff" F1 = "sd_diff"
putexcel A2 = "paired_ttest(a_hat_2-a_hat=0)" B2 = (`tstat') C2 = (`tpval') D2 = (`Ndiff') E2 = (`mdiff') F2 = (`sddiff')
putexcel A3 = "signrank(a_hat_2-a_hat=0)" B3 = (`sr_z') C3 = (`sr_p') D3 = (`Ndiff')

di as txt "Done. Output written to:"
di as txt "  ${OUTDIR}/Robust_2_compare.xlsx"
*******************************************************
