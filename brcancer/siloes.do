/*
siloes.do
IW 23jan2026
minor revision 9apr2026, 24apr2026
*/

// User-specific settings
cd C:\ian\git\TwoStagePoisson\brcancer
adopath ++ C:\ian\git\TwoStagePoisson\ado
set scheme mrc

version 18
cap log close
set linesize 100
log using siloes, text replace

version
which meta2p

local reps 50
set seed 461860

webuse brcancer, clear
stset rectime, fail(censrec)
stcox hormon, nohr
global b = _b[hormon]
global se = _se[hormon]

cap frame drop results
frame create results siloes rep b se
foreach siloes in 2 5 10 25 {
	frame post results (`siloes') (0) ($b) ($se)
	forvalues rep=1/`reps' {
		if `rep'==1 _dots 0, title("`siloes' siloes: simulation running (`reps' repetitions)")
		_dots `rep' 0
		cap drop random silo
		gen random = runiform()
		sort random
		gen silo = int(`siloes'*(_n-1)/_N)+1
		qui ipdmetan, study(silo) nograph: stcox hormon
		frame post results (`siloes') (`rep') (r(eff)) (r(se_eff))
		if `siloes'==25 & `rep'==11 qui save siloes`siloes'_`rep', replace
		* a fairly typical data set
	}
	frame results: label def siloes `siloes' "`siloes' siloes", modify	
}

frame results {
	sort siloes rep
	by siloes: ci mean b se
	label val siloes siloes
	save siloes, replace
}

frame results : ///
	scatter se b if rep==0, ms(X) msize(*2) ///
	|| scatter se b if rep>0, /*mlab(rep)*/ ms(oh) ///
	|| scatter se b if rep==11 & siloes==25, /*mlab(rep)*/ ms(O) msize(*1.5) ///
	by(siloes, note("")) legend(order(1 "Without splitting" 2 "2-stage Normal approx")) ///
	xline($b) yline($se) ///
	xtitle(Log hazard ratio) ytitle(Standard error) ///
	name(siloes, replace) $PPT


// explore a typical data set

use siloes25_11, clear

cap frame drop res25_11
frame create res25_11 str4 method b se
// one-stage
stcox hormon, nohr
frame post res25_11 ("1S") (_b[hormon]) (_se[hormon])

tab silo hormon if _d
byvar silo, b(hormon) se(hormon) unique gen: stcox hormon
egen d1 = sum(_d * hormon), by(silo)
egen d0 = sum(_d * (1-hormon)), by(silo)
egen p1 = sum(hormon), by(silo)
egen p0 = sum((1-hormon)), by(silo)
keep if !mi(B)
keep silo B S d? p?
rename B b
rename S se
format b se %6.3f
l, noo clean

// two-stage Normal
metan b se
frame post res25_11 ("2SN") (r(eff)) (r(se_eff))
// two-stage Poisson
meta2p b se, d(d1 d0) py(p1 p0)
frame post res25_11 ("2SP") (r(eff)) (r(se_eff))
meta2p b se, d(d1 d0) py(p1 p0) wt
frame post res25_11 ("2SPW") (r(eff)) (r(se_eff))

frame res25_11: l, noo clean

log close
