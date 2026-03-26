/*
twoSPcombine.do
IW 26mar2026
combine results of the two runs (same setting, different methods) in data twoSPcombine.dta
*/

// Check the two data sets agree

do twoSPsimrun2
simcombine, settings(simrun2_settings)
drop nobslow nobsupp px bx beta cens
isid studies aratio gamma tau rep
sort studies aratio gamma tau rep
su b*
* b1SC b1SCinter b1SWei b1SWei_ML b1SWei_MLc
recast double studies aratio gamma tau, force
replace tau = round(tau,1E-3)
save z2, replace

do twoSPsimrun
simcombine, settings(simrun_settings)
drop nobslow nobsupp px bx beta cens
isid studies aratio gamma tau rep
sort studies aratio gamma tau rep
su b*
* b1SC b2SN b2SPU b2SPW b2SN_REML b2SN_ML b2SPU_ML b2SPW_ML b2SPU_MLc b2SPW_MLc
recast double studies aratio gamma tau, force
replace tau = round(tau,1E-3)
save z1, replace

cf studies-error1SC using z2


// Merge the two data sets
use z2, clear
drop *1SC
merge 1:1 studies aratio gamma tau rep using z1
assert _merge==3
drop _merge

label data "Combined results from twoSPsimrun and twoSPsimrun2"
save twoSPcombine, replace

erase z1.dta
erase z2.dta
