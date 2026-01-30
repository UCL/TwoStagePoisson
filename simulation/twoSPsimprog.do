/*
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Simulation definition program (see twoSPsimrun for the main program)
IW 10sep2025
v0.2 IW 17sep2025
19sep2025: 1SC now stratifies, not adjusts (faster)
1oct2025: two programs, for revised simrun
*/
cap program drop twoSPsimgen
cap program drop twoSPsimana

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

prog def twoSPsimana, rclass

* descriptive outputs
local outputs n d singlezeroes doublezeroes sdb sdbw
* analysis outputs
foreach method in 1SC 2SN 2SPU 2SPW 2SN_REML 2SN_ML 2SPU_ML 2SPW_ML 2SPU_MLc 2SPW_MLc {
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

* reduce to two-stage data
keep if id==1

* CE Method 2: 2-stage Normal
timer on 2
cap metan b se, nograph model(iv)
timer off 2
if _rc==1 exit 1
if _rc==0 {
	local method 2SN
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = 0
	local error`method' = 0
}

* CE method 3: 2-stage Poisson, unweighted
timer on 3
cap meta2p b se, d(d1 d0) py(p1 p0) study(study)
timer off 3
if _rc==1 exit 1
if _rc==0 {
	local method 2SPU
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = 0
	local error`method' = r(error)
}

* CE method 4: 2-stage Poisson, weighted
timer on 4
cap meta2p b se, d(d1 d0) py(p1 p0) study(study) wt
timer off 4
if _rc==1 exit 1
if _rc==0 {
	local method 2SPW
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = 0
	local error`method' = r(error)
}

* RE Method 5: 2-stage Normal, REML
timer on 5
cap metan b se, nograph model(reml)
timer off 5
if _rc==1 exit 1
if _rc==0 {
	local method 2SN_REML
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = 0
}

* RE Method 6: 2-stage Normal, ML
timer on 6
cap metan b se, nograph model(ml)
timer off 6
if _rc==1 exit 1
if _rc==0 {
	local method 2SN_ML
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = 0
}

* RE method 7: 2-stage Poisson, unweighted
timer on 7
cap meta2p b se, d(d1 d0) py(p1 p0) study(study) re 
timer off 7
if _rc==1 exit 1
if _rc==0 {
	local method 2SPU_ML
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = r(error)
}

* RE method 8: 2-stage Poisson, weighted
timer on 8
cap meta2p b se, d(d1 d0) py(p1 p0) study(study) re wt
timer off 8
if _rc==1 exit 1
if _rc==0 {
	local method 2SPW_ML
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = r(error)
}

* RE method 9: 2-stage Poisson, unweighted, centred
timer on 9
cap meta2p b se, d(d1 d0) py(p1 p0) study(study) re centre
timer off 9
if _rc==1 exit 1
if _rc==0 {
	local method 2SPU_MLc
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = r(error)
}

* RE method 10: 2-stage Poisson, weighted, centred
timer on 10
cap meta2p b se, d(d1 d0) py(p1 p0) study(study) re wt centre
timer off 10
if _rc==1 exit 1
if _rc==0 {
	local method 2SPW_MLc
	local b`method' = r(eff)
	local s`method' = r(se_eff)
	local tau`method' = sqrt(r(tausq))
	local error`method' = r(error)
}

foreach output of local outputs {
	if !mi("``output''") return scalar `output' = ``output''
	else return scalar `output' = .
}

end

/************** end of twoSPsimana ******************/

