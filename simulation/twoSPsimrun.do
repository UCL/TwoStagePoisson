/*
twoSPsimrun.do
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Main program
2000 reps of 5 studies takes ~0.5-1hs; of 20 studies, 2-3h
IW Sep 2025
1oct2025 matches updated simrun.ado
19dec2025 for paper v0.3
29jan2026 for paper v0.4
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


// DESCRIPTIVE
simcombine 
foreach var in d n singlezeroes doublezeroes {
	replace `var'=`var'/studies
}
foreach var in singlezeroes doublezeroes {
	replace `var'=100*`var'
}
summ d n singlezeroes doublezeroes 
collapse (mean) d single doublez, by(gamma tau aratio)
order gamma tau aratio
format d %6.0f
format *zeroes %6.1f
l, noo clean abb(20)


// NON-CONVERGENCE
simcombine
mvpatterns b2SPU_ML b2SPU_MLc b2SPW_ML b2SPW_MLc


// COMPARE UNWEIGHTED AND WEIGHTED METHODS
simcombine
drop error* tau?*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(2SPU 2SPW 2SPU_ML 2SPW_ML) true(beta)
siman ana
siman nes pctbias modelse ///
	if inlist(method,"2SPU":methodlabel,"2SPW":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingCE,replace)
siman nes pctbias modelse ///
	if inlist(method,"2SPU_ML":methodlabel,"2SPW_ML":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingRE,replace)
table gamma studies if rep<0 & _perf=="pctbias", stat(max s) stat(min s) nformat(%6.3f)
* pctbias: Monte Carlo error varies from 0.9 to 2.7 with 5 studies and from 0.4 to 1.3 with 20 studies
* modelse: Monte Carlo error < 0.003


// COMPARE UNCENTRED AND CENTRED METHODS
simcombine
drop error* tau?*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(2SN_REML 2SN_ML 2SPU_MLc 2SPU_ML) true(beta)
siman ana
siman nes pctbias cover, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(blue = red =) lpat(solid dash solid dash) ///
	name(centring,replace)
* pctbias: Monte Carlo error varies from 0.9 to 2.9 with 5 studies and from 0.4 to 1.3 with 20 studies
* cover: Monte Carlo error varies from 0.3 to 0.9%


// SCATTER PLOT OF DIFFERENCES FROM 1SC: HOM, CE METHODS
simcombine
drop error* tau1* tau2* s1* s2*
foreach method in 2SN 2SPU {
	gen e`method' = b`method' - b1SC
}
reshape long e, i(studies aratio gamma tau rep) j(method) string
label var e "Difference from 1-stage Cox"
label var b1SC "1-stage Cox"
egen emean = mean(e), by(aratio method)
gen zero=0
recast double tau
replace tau = round(tau,.1)
keep if studies==5 & tau==0
*scatter e b1SC, mcol(blue) || line emean zero b1SC, lcol(blue black) by(method aratio gamma, legend(off) col(3) note("") t2title(gamma=-6                  gamma=-5                  gamma=-4)) ytitle("      2-stage Poisson                   2-stage Normal" "     aratio=2        aratio=1        aratio=2        aratio=1") subtitle("")
scatter e b1SC, mcol(blue) || line emean zero b1SC, lcol(blue black) by(method gamma aratio, legend(off) col(6) note("") t2title("gamma=-6                                      gamma=-5                                      gamma=-4" "aratio=1                aratio=2                aratio=1                aratio=2                aratio=1                aratio=2")) ytitle("    log HR, 2-stage Poisson      log HR, 2-stage Normal") xsize(9) ysize(5) subtitle("") name(discrep_ref,replace) xtitle("log HR, 1-stage Cox")

replace method = "2-stage Normal" if method=="2SN"
replace method = "2-stage Poisson" if method=="2SPU"
label var sdb "Standard deviation of point estimates"
scatter e sdb if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sd,replace)
scatter e sdbw if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sdw,replace)

* compare discriminatory ability of sdb and sdbw
gen discrep = abs(e)>=0.05 if !mi(e)
roccomp discrep sdb sdbw, graph name(roccomp,replace) summary plot1opts(ms(i)) plot2(ms(i))


// ANALYSE ESTIMATES OF MU
simcombine
drop error* tau?*
* drop all 2SP weighted
drop b2SPW s2SPW b2SPW_ML s2SPW_ML b2SPW_MLc s2SPW_MLc
* drop all 2SP ML uncentred
drop b2SPU_ML s2SPU_ML
* naming doesn't need to include U and c
rename (*2SPU*) (*2SP*)
rename (*MLc) (*ML)
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(1SC 2SN 2SP 2SN_REML 2SP_ML) true(beta)
siman ana

* bias & power, homogeneous & CE methods
siman tab pctbias empse cover power if tau==0 & method<=3, col(method gamma studies) row(_perf aratio) nformat(%6.1f)
siman nes pctbias empse cover power if tau==0 & method<=3, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red =) lpat(solid = = dash) name(homCE,replace)
* summarise MCSEs
table gamma studies if rep<0 & _perf=="pctbias" & tau==0 & method<=3, stat(max s) stat(min s) nformat(%6.1f)
* pctbias: MCSE varies from 0.9 to 2.7 with 5 studies and from 0.4 to 1.2 with 20 studies
* Empse: MCSE < 0.006
* Coverage: MCSE = 0.3 to 0.5%
* power: MCSE <= 1.1%

* bias, heterogeneous & RE methods
siman tab pctbias empse cover power if tau>0, col(method gamma studies) row(_perf aratio) nformat(%6.1f)
siman nes pctbias empse cover power if tau>0, ///
	legend(pos(3) col(1)) stagger(.08) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red blue red) lpat(solid = = dash =) name(het,replace)
siman nes power if tau>0 & inlist(method,4,5), ///
	legend(pos(3) col(1)) stagger(.08) ylab(0 50 100) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(blue red) lpat(dash =) name(het,replace)
* summarise MCSEs
table gamma studies if rep<0 & _perf=="pctbias" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Pctbias: Monte Carlo error varies from 1.4 to 2.9 with 5 studies and from 0.7 to 1.3 with 20 studies
table gamma studies if rep<0 & _perf=="empse" & tau>0, stat(min s) stat(max s) nformat(%6.3f)
* Empse: Monte Carlo error < 0.006
table gamma studies if rep<0 & _perf=="cover" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Coverage: Monte Carlo error = 0.4 to 1.0%
table gamma studies if rep<0 & _perf=="power" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Power: Monte Carlo error <= 1.1%


// ANALYSE ESTIMATES OF TAU
simcombine
drop error* b1* b2* s1* s2*
rename tau truetau
* drop all 2SP ML weighted or uncentred
drop tau2SPW tau2SPW_ML tau2SPW_MLc tau2SPU_ML 
* naming doesn't need to include U and c
rename (*2SPU*) (*2SP*)
rename (*MLc) (*ML)
siman setup, dgm(aratio gamma studies truetau) rep(rep) est(tau) method(2SN_REML 2SP_ML) true(truetau)
siman ana, max(100)
siman tab mean if truetau>0, row(truetau gamma aratio studies)
siman nes mean if truetau>0, ///
	legend(pos(3) col(1)) stagger(.05) ///
	$PPT dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(blue red) lpat(dash dash) name(tauhet,replace)
table gamma studies if rep<0 & _perf=="mean" & truetau>0, stat(max _se) stat(min _se) nformat(%6.3f)
* Mean: Monte Carlo error < 0.007

/* before using simrecreate,need:
recast double px bx beta
foreach var of varlist px bx beta {
	replace `var' = round(`var',1E-7)
}
*/