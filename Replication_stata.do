
********************************************************************************
* Replication.do  (Python -> Stata translation of /mnt/data/Replication.py)
* Notes:
* - Set your data folder below (two .dta files are required).
* - This script follows the numbered "Questions" in the Python notebook.
* - Graph commands are included where straightforward; feel free to comment out.
********************************************************************************

version 16.0
clear all
set more off
set linesize 255

********************************************************************************
* 0) PATHS (EDIT THESE)
********************************************************************************
* Folder containing:
*   - country_cleaned.dta
*   - country_sector_cleaned.dta
cd "/Users/raph/Library/Mobile Documents/com~apple~CloudDocs/M2/Panel Data/Project/Data replication/Data"
global DATA_DIR "`c(pwd)'"
global COUNTRY_DTA "$DATA_DIR/country_cleaned.dta"
global SECTOR_DTA  "$DATA_DIR/country_sector_cleaned.dta"

********************************************************************************
* Q1) Panel dimensions (country panel + country-sector panel)
********************************************************************************
di as txt "=============================="
di as txt "Question 1"
di as txt "=============================="

* Country panel: #countries + time span
use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

* Count unique countries
bys country: gen byte __tag_country = (_n==1)
quietly count if __tag_country
local n_countries = r(N)
drop __tag_country

quietly summarize year, meanonly
local t_min = r(min)
local t_max = r(max)
local max_t_span = `t_max' - `t_min' + 1

* Sector panel: #sectors + #country-sector individuals + consecutive obs
use "$SECTOR_DTA", clear

* Count unique sectors
bys industry: gen byte __tag_industry = (_n==1)
quietly count if __tag_industry
local n_sectors = r(N)
drop __tag_industry

* Create a country-sector id (numeric is safer than string)
egen long country_sector_id = group(country industry), label

* #Individuals (country-sector pairs)
bys country_sector_id: gen byte __tag_cs = (_n==1)
quietly count if __tag_cs
local n_individuals = r(N)

* Individuals with only 1 observation
bys country_sector_id: gen long __cs_N = _N
quietly count if __tag_cs & __cs_N==1
local n_single_obs = r(N)

* Individuals with at least 2 consecutive observations (year difference == 1 somewhere)
bys country_sector_id (year): gen int __year_diff = year - year[_n-1]
bys country_sector_id: egen byte __has_consec = max(__year_diff==1)
quietly count if __tag_cs & __has_consec==1
local consecutive_obs = r(N)

di as res "Number of Countries:                 `n_countries'"
di as res "Number of Sectors:                   `n_sectors'"
di as res "Total Country-Sector individuals:    `n_individuals'"
di as res "Period:                              `t_min' to `t_max' (Max T Span = `max_t_span')"
di as res "Individuals with only 1 observation: `n_single_obs'"
di as res "Individuals with ≥2 consecutive obs: `consecutive_obs'"

drop __tag_cs __cs_N __year_diff __has_consec

********************************************************************************
* Q4) Variable transformation & variance decomposition (country panel)
*     - pooled variance
*     - between variance (variance of country means)
*     - within variance (variance of deviations from country mean)
********************************************************************************
di as txt "=============================="
di as txt "Question 4"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

* Identify numeric variables (exclude year)
ds, has(type numeric)
local numvars `r(varlist)'
local numvars : list numvars - year

tempfile q4_results
postfile q4h str32 varname double pooled_var between_var within_var share_between share_within using `q4_results', replace

foreach v of local numvars {
    quietly count if !missing(`v')
    if r(N) < 2 continue

    * pooled variance
    quietly summarize `v' if !missing(`v'), detail
    local pooled = r(Var)
    if (`pooled'==0) continue

    * between (country means)
    tempvar mean_i within_i
    bys country: egen double `mean_i' = mean(`v')
    quietly summarize `mean_i' if !missing(`mean_i'), detail
    local between = r(Var)

    * within (de-meaned)
    gen double `within_i' = `v' - `mean_i'
    quietly summarize `within_i' if !missing(`within_i'), detail
    local within = r(Var)

    local share_b = 100*`between'/`pooled'
    local share_w = 100*`within'/`pooled'

    post q4h ("`v'") (`pooled') (`between') (`within') (`share_b') (`share_w')

    * keep transformed variables (as in Python block)
    cap drop `v'_between `v'_within
    gen double `v'_between = `mean_i'
    gen double `v'_within  = `within_i'

    drop `within_i'
}
postclose q4h

preserve
use `q4_results', clear
order varname pooled_var between_var within_var share_between share_within
format pooled_var between_var within_var %12.6g
format share_between share_within %9.2f
list, noobs abbreviate(24)
export delimited using "Question4_Variance_Decomposition.csv", replace
restore

********************************************************************************
* Q5) Distribution analysis (histogram + normal + KDE-ish) for:
*     dep_var = cotwo_total_per_cap
*     exp_var = log_gdp_per_cap
********************************************************************************
di as txt "=============================="
di as txt "Question 5"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

local dep_var cotwo_total_per_cap
local exp_var log_gdp_per_cap

foreach v in `dep_var' `exp_var' {
    di as txt "---- Summary (detail) for `v' ----"
    quietly summarize `v', detail
    di as res "N=" r(N) " mean=" %9.4f r(mean) " sd=" %9.4f r(sd) " skew=" %9.4f r(skewness) " kurt=" %9.4f r(kurtosis)

    * Histogram with normal overlay
    histogram `v', normal name(h_`v', replace)
    graph export "Question5_hist_`v'.png", replace

    * Kernel density (Stata's kdensity)
    kdensity `v', name(k_`v', replace)
    graph export "Question5_kdensity_`v'.png", replace
}

********************************************************************************
* Q6) First differences for all numeric vars (country panel)
********************************************************************************
di as txt "=============================="
di as txt "Question 6"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)}
else {
    gen long country_id = country
}
xtset country_id year

ds, has(type numeric)
local numvars `r(varlist)'
local numvars : list numvars - year

foreach v of local numvars {
    cap drop `v'_fd
    gen double `v'_fd = D.`v'
}

********************************************************************************
* Q7) Transformations for dep_var=cotwo_total_per_gdp and exp_var=fin_str2
*     within + first-difference + between (country means)
********************************************************************************
di as txt "=============================="
di as txt "Question 7"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}
xtset country_id year

local dep_var cotwo_total_per_gdp
local exp_var fin_str2

bys country_id: egen double dep_between = mean(`dep_var')   // line 252 (edited)
bys country_id: egen double exp_between = mean(`exp_var')   // line 253 (edited)

cap drop dep_within exp_within dep_fd exp_fd
gen double dep_within = `dep_var' - dep_between
gen double exp_within = `exp_var' - exp_between

sort country_id year   // <-- INSERT HERE (between old lines 257 and 258)
gen double dep_fd     = D.`dep_var'
gen double exp_fd     = D.`exp_var'



* Correlations (FD, Within, Between)
quietly corr exp_fd dep_fd if !missing(exp_fd, dep_fd)
matrix C = r(C)
local corr_fd = C[1,2]

quietly corr exp_within dep_within if !missing(exp_within, dep_within)
matrix Cw = r(C)
local corr_within = Cw[1,2]

preserve
collapse (mean) `exp_var' `dep_var', by(country)
quietly corr `exp_var' `dep_var'
matrix Cb = r(C)
local corr_between = Cb[1,2]
restore

di as res "First-Difference correlation: " %9.4f `corr_fd'
di as res "Within (one-way FE) corr:     " %9.4f `corr_within'
di as res "Between-country corr:         " %9.4f `corr_between'

********************************************************************************
* Q8) Joint relationship in FD space: scatter + fitted line
********************************************************************************
di as txt "=============================="
di as txt "Question 8"
di as txt "=============================="

twoway (scatter dep_fd exp_fd, msize(tiny)) ///
       (lfit dep_fd exp_fd), ///
       name(q8_scatter, replace) ///
       title("FD: dep_fd vs exp_fd") ///
       ytitle("Δ CO2/GDP") xtitle("Δ fin_str2")
graph export "Question8_FD_Correlation.png", replace

********************************************************************************
* Q9) Balanced panel selection + TWFE component for dep_var=cotwo_total_per_gdp
********************************************************************************
di as txt "=============================="
di as txt "Question 9"
di as txt "=============================="

use "$COUNTRY_DTA", clear

* Ensure numeric panel id (country_id)
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

* ------------------------------------------------------------
* 1) Find years with the maximum number of unique countries
*    (no egen nvals(); compute unique country count per year)
* ------------------------------------------------------------
sort year country_id
bys year country_id: gen byte __tag_year_country = (_n==1)
bys year: egen long __n_country_year = total(__tag_year_country)
drop __tag_year_country

quietly summarize __n_country_year, meanonly
local maxN = r(max)

* Keep only those "best-covered" years (max #countries)
keep if __n_country_year == `maxN'
drop __n_country_year

* List the balanced years and their count
levelsof year, local(bal_years)
local Tbal : word count `bal_years'
di as res "Balanced panel years: `bal_years'"
di as res "T (number of years): `Tbal'"

* ------------------------------------------------------------
* 2) Keep only countries that appear in ALL of those years
* ------------------------------------------------------------
bys country_id: gen long __ny = _N
keep if __ny == `Tbal'
drop __ny

* Report #countries in the balanced sample
bys country_id: gen byte __tag_c = (_n==1)
quietly count if __tag_c
local Nbal = r(N)
drop __tag_c
di as res "Balanced panel countries: `Nbal'"

* Set panel structure
sort country_id year
xtset country_id year

* ------------------------------------------------------------
* 3) TWFE transformation for dep_var
* ------------------------------------------------------------
local dep_var cotwo_total_per_gdp

egen double mu_grand = mean(`dep_var')
bys country_id: egen double mu_i = mean(`dep_var')
bys year:       egen double mu_t = mean(`dep_var')

capture drop dep_twfe
gen double dep_twfe = `dep_var' - mu_i - mu_t + mu_grand

* ------------------------------------------------------------
* 4) Time component: -mean_t + grand mean (plot)
* ------------------------------------------------------------
preserve
collapse (mean) mu_t mu_grand, by(year)
gen double time_component = -mu_t + mu_grand
twoway line time_component year, ///
    name(q9_time, replace) ///
    title("Time component: -mean_t + grand mean") ///
    ytitle("Value") xtitle("Year")
graph export "Question9_Time_Component.png", replace
restore


********************************************************************************
* Q10) Descriptive stats in balanced panel
********************************************************************************
di as txt "=============================="
di as txt "Question 10"
di as txt "=============================="

* (we are still in the balanced sample from Q9)
tabstat cotwo_total_per_gdp fin_str2 log_gdp_per_cap dep_twfe, ///
    stat(n mean sd min p25 p50 p75 max) columns(statistics)

********************************************************************************
* Q11) TWFE for y=cotwo_total_per_gdp and x=log_gdp_per_cap, then by-country rho + slope
********************************************************************************
di as txt "=============================="
di as txt "Question 11"
di as txt "=============================="

local yvar cotwo_total_per_gdp
local xvar log_gdp_per_cap

* TWFE for both
egen double y_mu = mean(`yvar')
bys country: egen double y_mu_i = mean(`yvar')
bys year:    egen double y_mu_t = mean(`yvar')
cap drop y_twfe
gen double y_twfe = `yvar' - y_mu_i - y_mu_t + y_mu

egen double x_mu = mean(`xvar')
bys country: egen double x_mu_i = mean(`xvar')
bys year:    egen double x_mu_t = mean(`xvar')
cap drop x_twfe
gen double x_twfe = `xvar' - x_mu_i - x_mu_t + x_mu

tempfile q11
postfile q11h str32 country double rho sd_x sd_y reg_coef using `q11', replace

* Loop by numeric country_id to avoid string/numeric type-mismatch issues
levelsof country_id, local(cid_list)
foreach cid of local cid_list {
    preserve
    keep if country_id==`cid'
    * derive a readable country label/name
    local cname "`cid'"
    capture confirm string variable country
    if !_rc {
        quietly levelsof country if country_id==`cid', local(_nm)
        local cname : word 1 of `_nm'
    }
    else {
        capture local cname : label (country_id) `cid'
    }

    quietly corr x_twfe y_twfe if !missing(x_twfe, y_twfe)
    matrix CC = r(C)
    local rho = CC[1,2]
    quietly summarize x_twfe, meanonly
    local sdx = r(sd)
    quietly summarize y_twfe, meanonly
    local sdy = r(sd)
    local slope = cond(`sdy'!=0, `rho'*(`sdx'/`sdy'), 0)
    post q11h ("`cname'") (`rho') (`sdx') (`sdy') (`slope')
    restore
}
postclose q11h

preserve
use `q11', clear
gsort -rho
format rho sd_x sd_y reg_coef %9.4f
list, noobs abbreviate(24)
export delimited using "Question11_TWFE_Correlations.csv", replace
restore

********************************************************************************
* Q12) Unbalanced TWFE residuals via FE regression with year dummies:
*     residuals from: var = country FE + year FE + e  (works for unbalanced panels)
********************************************************************************
di as txt "=============================="
di as txt "Question 12"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}
bys country: gen long __N = _N
keep if __N>1
drop __N
xtset country_id year

* y and x
local y_col cotwo_total_per_gdp
local x_col log_gdp_per_cap

xtreg `y_col' i.year, fe
predict double y_twfe_unbal, e

xtreg `x_col' i.year, fe
predict double x_twfe_unbal, e

summ y_twfe_unbal x_twfe_unbal

********************************************************************************
* Q14–Q15) (Core parts) Build transformations for y=cotwo_total_per_gdp and x=fin_str2
*         - between, within, twfe, fd
*         - compute univariate stats with standardized extremes (min/max z)
********************************************************************************
di as txt "=============================="
di as txt "Questions 14–15 (core)"
di as txt "=============================="

use "$COUNTRY_DTA", clear

* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

* Keep only countries with >1 observation
sort country_id year
bys country_id: gen long __N = _N
keep if __N>1
drop __N

* Declare panel + enforce sort for time-series operators
sort country_id year
xtset country_id year

local y_raw cotwo_total_per_gdp
local x_raw fin_str2

* -----------------------------
* Between + Within
* -----------------------------
sort country_id year
bys country_id: egen double y_between = mean(`y_raw')
bys country_id: egen double x_between = mean(`x_raw')

gen double y_within = `y_raw' - y_between
gen double x_within = `x_raw' - x_between

* -----------------------------
* TWFE residuals (unbalanced ok): country FE + year FE
* -----------------------------
xtreg `y_raw' i.year, fe
predict double y_twfe, e

xtreg `x_raw' i.year, fe
predict double x_twfe, e

* -----------------------------
* First differences (THIS is where you were getting "not sorted")
* -----------------------------
sort country_id year
gen double y_fd = D.`y_raw'
gen double x_fd = D.`x_raw'

* -----------------------------
* Univariate stats table (Mean, SD, quartiles, z-min, z-max, mean-median)
* -----------------------------
tempfile q15
postfile q15h str10 var str10 trans double mean sd zmin q1 med q3 zmax mean_minus_med using `q15', replace

foreach var in y x {
    foreach tr in between within twfe fd {
        local vname `var'_`tr'
        quietly summarize ``vname'', detail
        local mu = r(mean)
        local s  = r(sd)
        local mn = r(min)
        local mx = r(max)
        local p25 = r(p25)
        local p50 = r(p50)
        local p75 = r(p75)
        local zmin = cond(`s'!=0, (`mn'-`mu')/`s', .)
        local zmax = cond(`s'!=0, (`mx'-`mu')/`s', .)
        local mdiff = `mu' - `p50'
        post q15h ("`var'") ("`tr'") (`mu') (`s') (`zmin') (`p25') (`p50') (`p75') (`zmax') (`mdiff')
    }
}
postclose q15h

use `q15', clear
format mean sd zmin zmax mean_minus_med %9.4f
list, noobs abbreviate(24)
export delimited using "Question15_Descriptive_Stats.csv", replace


********************************************************************************
* Q20) Model comparison (Between, FE, Mundlak-RE, TWFE, FD)
********************************************************************************
di as txt "=============================="
di as txt "Question 20"
di as txt "=============================="

use "$COUNTRY_DTA", clear

* Ensure numeric panel id (country_id)
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}

* Keep only countries with >1 observation
sort country_id year
bys country_id: gen long __N = _N
keep if __N>1
drop __N

xtset country_id year
sort country_id year

local dep cotwo_total_per_gdp
local X "fin_str2 log_gdp_per_cap credit stock pop_mil"

* Comparable sample (drop missing in all vars)
keep if !missing(`dep', fin_str2, log_gdp_per_cap, credit, stock, pop_mil)

* ------------------------------------------------------------
* Between (BE): do it manually so we can use robust/cluster SEs
* ------------------------------------------------------------
preserve
collapse (mean) `dep' fin_str2 log_gdp_per_cap credit stock pop_mil, by(country_id)

* Robust SE:
reg `dep' fin_str2 log_gdp_per_cap credit stock pop_mil, vce(robust)
estimates store BE

* (If you prefer cluster by country_id, use this instead:)
* reg `dep' fin_str2 log_gdp_per_cap credit stock pop_mil, vce(cluster country_id)
* estimates store BE

restore

* One-way FE (country FE)
xtreg `dep' `X', fe vce(cluster country_id)
estimates store FE1

* Mundlak (RE with country means of X)
foreach v in fin_str2 log_gdp_per_cap credit stock pop_mil {
    bys country_id: egen double mean_`v' = mean(`v')
}
xtreg `dep' `X' mean_fin_str2 mean_log_gdp_per_cap mean_credit mean_stock mean_pop_mil, re vce(cluster country_id)
estimates store MUNDLAK

* TWFE (country FE + year FE)
xtreg `dep' `X' i.year, fe vce(cluster country_id)
estimates store TWFE

* First difference (FD)
xtreg `dep' `X' i.year, fd vce(cluster country_id)
estimates store FD

* If you have esttab installed:
* esttab BE FE1 MUNDLAK TWFE FD, se ar2


********************************************************************************
* Q21) Anderson–Hsiao ARDL(1,1) in first differences + IV
********************************************************************************
di as txt "=============================="
di as txt "Question 21 (Anderson–Hsiao)"
di as txt "=============================="

use "$COUNTRY_DTA", clear
* Ensure numeric panel id (country_id) for xt* commands
capture drop country_id
capture confirm numeric variable country
if _rc {
    encode country, gen(country_id)
}
else {
    gen long country_id = country
}
bys country: gen long __N = _N
keep if __N>1
drop __N
xtset country_id year

local y cotwo_total_per_gdp
local x fin_str2
local controls "log_gdp_per_cap credit stock pop_mil"

* OLS ARDL in differences (with year FE)
reg D.`y' L.D.`y' D.`x' L.D.`x' D.(log_gdp_per_cap credit stock pop_mil) i.year, vce(cluster country_id)
estimates store ARDL_OLS

* IV (endogenous: L.D.y and D.x; instruments: L2.y and L2.x)
ivregress 2sls D.`y' L.D.`x' D.(log_gdp_per_cap credit stock pop_mil) i.year ///
    (L.D.`y' D.`x' = L2.`y' L2.`x'), vce(cluster country_id)
estimates store ARDL_IV

* First-stage and tests (built-in)
estat firststage
estat endogenous
capture noisily estat overid

* Q21.7 IRF for t=1..4 based on IV coefficients
scalar beta_y = _b[L.D.`y']
scalar beta_1 = _b[D.`x']
scalar beta_2 = _b[L.D.`x']

clear
set obs 4
gen int t = _n
gen double irf = .
replace irf = beta_1                                   if t==1
replace irf = beta_y*beta_1 + beta_2                   if t==2
replace irf = (beta_y^2)*beta_1 + beta_y*beta_2        if t==3
replace irf = (beta_y^3)*beta_1 + (beta_y^2)*beta_2    if t==4

twoway line irf t, ///
    title("IRF: Δy response to 1-unit shock in Δx") ///
    ytitle("Impact on Δy") xtitle("Periods (years)") ///
    name(q21_irf, replace)
graph export "Question21_7_IRF_Plot.png", replace

* Q21.8 Long-run propensity: (beta1 + beta2) / (1 - beta_y)
scalar beta_LT = (beta_1 + beta_2) / (1 - beta_y)
di as txt "----------------------------------------"
di as txt "Question 21.8 Long-Run Propensity (LRP)"
di as txt "beta_1 = " %9.6f beta_1
di as txt "beta_2 = " %9.6f beta_2
di as txt "beta_y = " %9.6f beta_y
di as txt "beta_LT= " %9.6f beta_LT
di as txt "----------------------------------------"

********************************************************************************
* End
********************************************************************************
