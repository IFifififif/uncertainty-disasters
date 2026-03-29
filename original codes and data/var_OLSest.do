**************************************************************
* var_OLSest.do
*
* This STATA code estimates the no country FE coefficient matrices
* for the disaster event restrictions VAR in 
* "Using Disasters to Estimate the Impact of Uncertainty"
* by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry 
*
* Initially run on STATA/MP 17.0
*
* Questions to Stephen Terry
* stephenjamesterry@gmail.com
**************************************************************

local data "./"
local results "./NO_COUNTRY_FE/"

use "`data'Dates_and_Data.dta", clear

*number of VAR lags (quarters)
local nlags = 12
scalar nlags_sc = 12

*bookkeeping to set up quarter indexes
gen yqind = (year-1970)*4+(quarter-1)

*encode countries as integers
encode country, gen(cnum)

*create panel structure with country x quarter
xtset cnum yqind

*set up the 3 variables in the VAR
gen est_y = ydgdp
gen est_first = l.ret
gen logvol = log(vol)
gen est_second = l.logvol

summ est_first 
replace est_first = est_first/r(sd)

summ est_second
replace est_second = est_second/r(sd)


*create lagged variables for the VAR OLS estimation
foreach lagnum of numlist 1/`nlags' {

gen est_y_`lagnum' = l`lagnum'.est_y
gen est_first_`lagnum' = l`lagnum'.est_first
gen est_second_`lagnum' = l`lagnum'.est_second
	
	
}

*create the common VAR sample
gen VARsamp = (1==1)
replace VARsamp = VARsamp&!missing(est_y)
replace VARsamp = VARsamp&!missing(est_first)
replace VARsamp = VARsamp&!missing(est_second)

foreach lagnum of numlist 1/`nlags' {

replace VARsamp = VARsamp & !missing(est_y_`lagnum')
replace VARsamp = VARsamp & !missing(est_first_`lagnum')
replace VARsamp = VARsamp & !missing(est_second_`lagnum')

	
	
}

foreach var of varlist revshock polshock tershock natshock {
	gen `var'_jump = `var'*(jump-1)*(jump>0)
	replace `var'_jump = 0 if missing(jump)
	
	gen `var'_postpred = `var'*post_pred
	replace `var'_postpred = 0 if missing(post_pred)

	
}

reghdfe est_y est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_y)
replace VARsamp = VARsamp & !missing(res_y)
reghdfe est_first est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_first)
replace VARsamp = VARsamp & !missing(res_first)
reghdfe est_second est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_second)
replace VARsamp = VARsamp & !missing(res_second)
drop res_y res_first res_second

*estimate the VAR
reghdfe est_y est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_y)
scalar lnL_y = e(ll)
matrix Avec_y = e(b)
mat2txt,matrix(Avec_y) saving("`results'Avec_y.txt") replace

reghdfe est_first est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_first)
scalar lnL_first = e(ll)
matrix Avec_first = e(b)
mat2txt,matrix(Avec_first) saving("`results'Avec_first.txt") replace

reghdfe est_second est_y_* est_first_* est_second_* if VARsamp, absorb(yqind) res(res_second)
scalar lnL_second = e(ll)
matrix Avec_second = e(b)
mat2txt,matrix(Avec_second) saving("`results'Avec_second.txt") replace


keep if VARsamp

export delimited cnum yqind est_y est_first est_second res_y res_first res_second jump post_pred *shock *_jump *_postpred  using "`results'VARout.csv", nolabel replace
