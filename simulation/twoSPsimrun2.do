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

