/*
*! v0.5.1 IW 17feb2026
	options for Poisson must be included in new poissonoptions()
	this means getting an option wrong (e.g. random instead of re) gives a clearer error message
v0.5 IW 10dec2025
	add centre option: likely to outperform nocentre for REML
v0.4 IW 25sep2025
	wt fails with zeroes
v0.3 IW 19sep2025
	use pyears for single-zero studies
v0.2 IW 10sep2025
*/
prog def meta2p, rclass
* run 2-stage Poisson-approx MA (CE and RE)
* NB b is assumed to be a log HR or log RR for group 1 vs group 2
syntax varlist(min=2 max=2) [if] [in], d(varlist min=2 max=2) [Study(varname) PYears(varlist min=2 max=2) wt RE CENtre VERBose pause list irr eform POISSONoptions(string) ///
debug /// undocumented
]
local est : word 1 of `varlist'
local se : word 2 of `varlist'
local d1 : word 1 of `d'
local d0 : word 2 of `d'
if !mi("`verbose'") local noi noi
if !mi("`irr'") local eform eform
* END OF PARSING

preserve
marksample touse
qui keep if `touse'
if missing("`study'") {
	tempvar study
	gen `study'=_n
}

qui count if max(`d0',`d1')==0
local ndoublezero = r(N)
qui count if min(`d0',`d1')==0
local nsinglezero = r(N) - `ndoublezero'
if `ndoublezero' {
	di as text "Ignoring " as result `ndoublezero' as text " double-zero studies"
}
if `nsinglezero' {
	di as text "Found " as result `nsinglezero' as text " single-zero studies: comparing events/person-years " _c
	if !mi("`pyears'") di as text "using reported person-years"
	else di as text "assuming equal person-years"
}
if !missing("`pyears'") { // these are observed pyears
	local py1 : word 1 of `pyears'
	local py0 : word 2 of `pyears'
	local pyratio `py0'/`py1'
}
else local pyratio 1

// Create Poisson pseudo-person-years and reshape
tempvar exp d z // exp holds pseudo-pyears
gen `exp'1 = 1
qui gen `exp'0 = (`d0'/`d1')*exp(`est') if min(`d0',`d1')>0
if `nsinglezero' qui replace `exp'0 = `pyratio' if min(`d0',`d1')==0 
rename `d1' `d'1
rename `d0' `d'0
qui keep `study' `d'1 `d'0 `exp'1 `exp'0 `se' `pyears'
if !mi("`wt'") {
	tempvar wtvar
	qui gen `wtvar' = (1/`d'1+1/`d'0)/`se'^2
	qui replace `wtvar' = 1 if min(`d'0,`d'1)==0 
	local wtexp [iweight=`wtvar']
	cap assert !mi("`wtvar'")
	char `wtvar'[varname] "Weight"
}
drop `se'
qui reshape long `d' `exp', i(`study') j(pooled)
rename `d' Meta

// Centring
if !mi("`centre'") {
	qui replace pooled = pooled - 1/(1+`pyratio')
}

// Fit Poisson model
local poissoncmd poisson Meta pooled i.`study' `wtexp', exp(`exp')
if !mi("`re'") local poissoncmd me`poissoncmd' || `study':pooled, nocons
local poissoncmd `poissoncmd' `poissonoptions'
if !mi("`list'") {
	di as text "Data for Poisson model:" _c
	char `exp'[varname] "Pyears"
	char pooled[varname] "Group"
	char Meta[varname] "Events"
	l `study' pooled Meta `exp' `wtvar' , sepby(`study') subvarname
}
if !mi("`debug'") di as input "Debug: `poissoncmd'"
if !mi("`pause'") {
	global F9 `poissoncmd'
	pause
}
cap `noi' `poissoncmd'
if _rc local error 2
else if !e(converged) local error 1
else local error 0
if `error'==2 di as error "Warning: [me]poisson ended with error " _rc
if `error'==1 di as error "Warning: [me]poisson did not converge"
if `error' {
	di as error `"command was: `poissoncmd'"'
	exit _rc
}

// Display results
local model = cond(mi("`re'"),"common","random")
local weighting = cond(mi("`wt'"),"unweighted","weighted")
di as text "Two-stage `model'-effects meta-analysis by `weighting' Poisson approximation:" _c
if !mi("`re'") local nlcom2 (tausq:_b[/var(pooled[`study'])])
nlcom (pooled:_b[pooled]) `nlcom2', noheader `eform'

// Return results, using same scalar names as metan
mat b=r(b)
mat V=r(V)
return scalar eff = b[1,"pooled"]
return scalar se_eff = sqrt(V["pooled","pooled"])
if !mi("`re'") return scalar tausq = b[1,"tausq"]
return scalar error = `error'
end

