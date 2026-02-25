/*
twoSPsimrun.do
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Main program
2000 reps of 5 studies takes ~0.5-1hs; of 20 studies, 2-3h
IW Sep 2025
1oct2025 matches updated simrun.ado
19dec2025 for paper v0.3
29jan2026 for paper v0.4
25feb2026 moved graphs and tables to twoSPsimresults.do
*/

// User-specific settings
cd C:\ian\git\TwoStagePoisson\simulation
adopath ++ C:\ian\git\TwoStagePoisson\ado
set scheme mrc

// SET UP DGMS
prog drop _all
cap frame create simrun_settings
frame change simrun_settings
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
gsort studies -tau gamma aratio // start with k=5 (faster), and tau=0.3 (more problematic)

* load the simulation program
do twoSPsimprog

simsetup, name(twoSPsimprog) prog(twoSPsimgen twoSPsimana)
timer clear
simrun, append
timer list


