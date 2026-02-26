/*
twoSPsimrun2.do
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Main program
Part 2: twoSPsimprog2.do changes analysis only to include more one-stage models:
	Cox CE with interaction (1SCinter)
	Weibull CE (1SWei)
	Weibull RE (1SWei_ML)
IW 25feb2026
*/

// User-specific settings
cd C:\ian\git\TwoStagePoisson\simulation
adopath ++ C:\ian\git\TwoStagePoisson\ado
set scheme mrc

// SET UP DGMS
prog drop _all
cap frame create simrun2_settings
frame change simrun2_settings
clear 
makegrid tau, values(0 0.3)
makegrid gamma, values(-6 -5 -4)
makegrid aratio, values(1 2)
makegrid studies, values(5 20)
gen bx = .4 
gen beta = -0.3
gen cens = 5
gen nobslow = 300
gen nobsupp = 1000
gen px = .4

gen reps = 2000
gen seed = 42
gsort studies -tau -gamma aratio // start with k=5 (faster) and tau=0.3 (more problematic)

* load the simulation program
do twoSPsimprog2

simsetup, name(twoSPsim2) folder(simrun2_results) prog(twoSPsimgen twoSPsimana2)
timer clear
simrun, append settings(simrun2_settings)
timer list

// INITIAL RESULTS 
simcombine in 1/8, settings(simrun2_settings)
tab1 studies tau gamma aratio
bysort studies tau gamma aratio: su b1S* tau1SWei_ML*
ci mean b1S*
foreach method in Cinter Wei Wei_ML Wei_MLc {
	qui gen diff`method' = b1S`method' - b1SC
	label var diff`method' "b1S`method' - b1SC"
	scatter diff`method' b1SC, by(studies tau gamma aratio) name(diff`method', replace)
}

simcombine in 1/8, settings(simrun2_settings)
keep studies tau aratio gamma rep b1S* s1S*
siman setup, dgm(studies tau gamma aratio) method(1SC 1SCinter 1SWei 1SWei_ML 1SWei_MLc) rep(rep) est(b) se(s) true(-.3)
siman scatter
siman ana
siman table pctbias, row(tau gamma aratio) nformat(%6.0f)
siman nes pctbias, stagger(.1) legend(pos(3) col(1)) name(mu)

simcombine in 1/8, settings(simrun2_settings)
keep studies tau aratio gamma rep tau1SWei_ML*
rename tau tautrue
siman setup, dgm(studies tautrue gamma aratio) method(1SWei_ML 1SWei_MLc) rep(rep) est(tau) true(tautrue)
siman swarm
siman ana, max(20)
siman table bias, row(tautrue gamma aratio) nformat(%6.2f)
siman nes bias, stagger(.1) legend(pos(3) col(1)) name(tau, replace)

