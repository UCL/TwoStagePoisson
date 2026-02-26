/*
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Simulation definition program (see twoSPsimrun for the main program)
Part 2: change analysis only to include more one-stage models:
	Cox CE with interaction (1SCinter)
	Weibull CE (1SWei)
	Weibull RE (1SWei_ML)
	Weibull RE + centred trt (1SWei_MLc)
IW 25feb2026
*/
cap program drop twoSPsimgen
cap program drop twoSPsimana2

program define twoSPsimgen, rclass

// FOR SIMRUN
local inputs studies nobslow nobsupp px aratio gamma bx beta tau cens
if "`0'"=="?" {
	return local inputs `inputs'
	exit
}

// GENERATE THE DATA
version 14
syntax, studies(int) nobslow(int) nobsupp(int) px(real) aratio(real) gamma(real) bx(real) beta(real) tau(real) cens(real)
assert `studies'>1
assert inrange(`px',0,1)
assert `aratio'>0
clear
qui set obs `studies'
gen study = _n
gen nobs = runiformint(`nobslow',`nobsupp')
gen bz = rnormal(`beta',`tau')
qui expand nobs
drop nobs
sort study
by study: gen id = _n
gen x = runiform()<`px'
gen z = runiform()>1/(1+`aratio')
gen t = -log(runiform())/exp(`gamma'+`bx'*x+bz*z)
gen d = t<=`cens'
qui replace t = `cens' if !d
drop bz
stset t, failure(d)

* Form 2-stage data
byvar study, b(z) se(z) generate unique: stcox z x
rename Bz_ b
rename Sz_ se
assert mi(b)==(id>1)
egen d1 = sum(d*z), by(study)
egen d0 = sum(d*(1-z)), by(study)
egen p1 = sum(t*z), by(study)
egen p0 = sum(t*(1-z)), by(study)
foreach var of varlist d? p? {
	qui replace `var' = . if id>1
}
* handle single-zero studies
qui replace b = ln(d1+.5)-ln(p1)-ln(d0+.5)+ln(p0) if min(d0,d1)==0 
qui replace se = sqrt( 1/(d1+.5) + 1/(d0+.5) )    if min(d0,d1)==0
* handle double-zero studies
qui replace b = . if max(d0,d1)==0
qui replace se = . if max(d0,d1)==0

end

/************** end of twoSPsimgen ******************/

/************** start of twoSPsimana ******************/

prog def twoSPsimana2, rclass

* descriptive outputs
local outputs n d singlezeroes doublezeroes sdb sdbw
* analysis outputs
foreach method in 1SC 1SCinter 1SWei 1SWei_ML 1SWei_MLc {
	local outputs `outputs' b`method' s`method' tau`method' error`method' 
}
if "`0'"=="?" {
	return local outputs `outputs'
	exit
}

// DESCRIPTIVE STATISTICS
sum d, meanonly
local n=r(N)
local d=r(sum)
qui count if d0==0 & d1==0 & id==1
local doublezeroes = r(N)
qui count if (d0==0 | d1==0) & id==1
local singlezeroes = r(N)-`doublezeroes'
* var(b)
qui sum b
local sdb = r(sd)
qui sum b [w=1/se^2]
local sdbw = r(sd)

// ANALYSE THE DATA
* Method 1: one-stage Cox
timer on 1
cap stcox z x, strata(study)
timer off 1
if _rc==1 exit 1
if _rc==0 {
	local method 1SC
	local b`method' = _b[z]
	local s`method' = _se[z]
	local tau`method' = 0
	local error`method' = !e(converged)
}

* Method 11: one-stage Cox with interaction
timer on 11
cap stcox z i.study##x, strata(study)
timer off 11
if _rc==1 exit 1
if _rc==0 {
	local method 1SCinter
	local b`method' = _b[z]
	local s`method' = _se[z]
	local tau`method' = 0
	local error`method' = !e(converged)
}

* Method 12: one-stage Weibull
timer on 12
cap streg z x, strata(study) dist(weibull)
timer off 12
if _rc==1 exit 1
if _rc==0 {
	local method 1SWei
	local b`method' = _b[z]
	local s`method' = _se[z]
	local tau`method' = 0
	local error`method' = !e(converged)
}

* Method 13: one-stage Weibull with RE
timer on 13
cap mestreg z x i.study || study:z, nocons dist(weibull)
timer off 13
if _rc==1 exit 1
if _rc==0 {
	local method 1SWei_ML
	local b`method' = _b[z]
	local s`method' = _se[z]
	local tau`method' = sqrt(_b[/var(z[study])])
	local error`method' = !e(converged)
}

* Method 14: one-stage Weibull with RE, treatment centred
timer on 14
summ z, meanonly
gen cz = z - r(mean)
cap mestreg cz x i.study || study:cz, nocons dist(weibull)
timer off 14
if _rc==1 exit 1
if _rc==0 {
	local method 1SWei_MLc
	local b`method' = _b[cz]
	local s`method' = _se[cz]
	local tau`method' = sqrt(_b[/var(cz[study])])
	local error`method' = !e(converged)
}

foreach output of local outputs {
	if !mi("``output''") return scalar `output' = ``output''
	else return scalar `output' = .
}

end

/************** end of twoSPsimana ******************/

