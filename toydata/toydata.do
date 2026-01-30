/*
"toydata.do"
IW 5/9/2025
for github 30/1/2026
*/

// User-specific settings
cd C:\ian\git\TwoStagePoisson\toydata
adopath ++ C:\ian\git\TwoStagePoisson\ado
set scheme mrc

// Compute PLLFs: a bit slow
prog drop _all
use toydata, clear
forvalues i=1/2 {
	poisson d z if study==`i', exposure(py) 
	local b`i'=_b[z]
	local se`i'=_se[z]
	pllf, profile(z) range(-3 2) nograph diff norm n(500) ///
		gen(beta`i' pll`i' norm`i'): ///
		poisson d z if study==`i', exposure(py) 
	label var pll`i' "Study `i', exact PLLF"
	label var norm`i' "Study `i', normal approx"
	label var beta`i' "beta: log rate ratio for treatment"
}
assert beta1==beta2
drop study z d py beta2
rename beta1 beta
save toydata_pllf, replace

// Rearrange data
use toydata_pllf, clear
gen pll12=pll1+pll2
gen norm12=norm1+norm2
summ pll12
replace pll12=pll12-r(max)
summ norm12
replace norm12=norm12-r(max)

reshape long pll norm, i(beta) j(method)
label var pll "Exact PLLF"
label var norm "Normal approx"
label def method 1 "Study 1" 2 "Study 2" 12 "Pooled"
label val method method

// Draw graph
foreach var in pll norm {
	replace `var'=. if `var'<-15
}
keep if inrange(beta,-3,1)
line pll norm beta, lc(black =) lp(solid dash) lw(*2 =) ytitle(Relative profile log-likelihood,size(small)) legend(row(1) size(small)) xsize(6) ysize(7) xtitle(,size(small)) ytitle(,size(small)) ylab(,labsize(small)) xlab(,labsize(small)) name(toydataPLL,replace) by(method, note("") col(1) legend(pos(6))) subtitle(,size(small)) xline(-0.355, lcol(gray)) saving(toydata_pllf,replace)
