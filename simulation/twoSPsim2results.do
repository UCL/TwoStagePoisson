/*
twoSPsimrun2: initial results 
twoSPsim2results.do
*/
simcombine, settings(simrun2_settings)
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

