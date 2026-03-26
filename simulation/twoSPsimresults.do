/*
twoSPsimresults.do
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Main results
*/

// DESCRIPTIVE
use twoSPcombine, clear 
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
use twoSPcombine, clear
mvpatterns b1*
mvpatterns b2*


// COMPARE 1SC METHODS
use twoSPcombine, clear
gen absdiff=abs(b1SWei - b1SCox)
su absdiff
* 1SWei and 1SC never differ by more than 0.004

drop error* tau?*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(1SCox 1SCoxi 1SWei) true(beta)
siman ana, max(15)
siman nes pctbias, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black green purple) lpat(solid = =) ///
	name(compare1SC,replace) yla(-15(5)5)
siman nes modelse, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black green purple) lpat(solid = =) ///
	name(compare1SC,replace) yla(0 .1 .2 .3)
* the 3 methods have indistinguishable bias and modelSE
table studies method if rep==-4, stat(min s) stat(max s)
table studies method if rep==-6, stat(min s) stat(max s)


// COMPARE UNWEIGHTED AND WEIGHTED 2SP METHODS
use twoSPcombine, clear
drop error* tau?*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(2SPU 2SPW 2SPU_REc 2SPW_REc) true(beta)
siman ana
siman nes pctbias ///
	if inlist(method,"2SPU":methodlabel,"2SPW":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingCE,replace) yla(-15(5)5)
siman nes modelse ///
	if inlist(method,"2SPU":methodlabel,"2SPW":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingCE,replace) yla(.1 .2 .3)
siman nes pctbias ///
	if inlist(method,"2SPU_REc":methodlabel,"2SPW_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingRE,replace) yla(-15(5)5)
siman nes modelse ///
	if inlist(method,"2SPU_REc":methodlabel,"2SPW_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingRE,replace) yla(.1 .2 .3 .4)
table gamma studies if rep<0 & _perf=="pctbias", stat(max s) stat(min s) nformat(%6.3f)
* pctbias: Monte Carlo error varies from 0.9 to 2.7 with 5 studies and from 0.4 to 1.3 with 20 studies
* modelse: Monte Carlo error < 0.003


// COMPARE UNCENTRED AND CENTRED METHODS
use twoSPcombine, clear
drop error* tau?*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(1SWei_REc 1SWei_REu 2SPU_REc 2SPU_REu) true(beta)
siman ana
siman nes pctbias ///
	if inlist(method,"2SPU_REu":methodlabel,"2SPU_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(centring2SPU,replace) yla(-15(5)5)
siman nes cover ///
	if inlist(method,"2SPU_REu":methodlabel,"2SPU_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(centring2SPU,replace) yla(80(5)100)
siman nes pctbias ///
	if inlist(method,"1SWei_REu":methodlabel,"1SWei_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black =) lpat(solid dash) ///
	name(centring1SWei,replace) yla(-15(5)5)
siman nes cover ///
	if inlist(method,"1SWei_REu":methodlabel,"1SWei_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black =) lpat(solid dash) ///
	name(centring1SWei,replace) yla(80(5)100)
table studies method if rep==-4, stat(min s) stat(max s)
* pctbias: Monte Carlo error varies from 0.9 to 2.9 with 5 studies and from 0.4 to 1.3 with 20 studies
table studies method if rep==-13, stat(min s) stat(max s)
* cover: Monte Carlo error varies from 0.4 to 0.9%.


// SCATTER PLOT OF DIFFERENCES FROM 1SC: HOM, CE METHODS
use twoSPcombine, clear
drop error* tau1* tau2* s1* s2*
foreach method in 2SN 2SPU {
	gen e`method' = b`method' - b1SCox
}
reshape long e, i(studies aratio gamma tau rep) j(method) string
label var e "Difference from 1-stage Cox"
label var b1SCox "1-stage Cox"
egen emean = mean(e), by(aratio method)
gen zero=0
recast double tau
replace tau = round(tau,.1)
keep if studies==5 & tau==0
*scatter e b1SC, mcol(blue) || line emean zero b1SC, lcol(blue black) by(method aratio gamma, legend(off) col(3) note("") t2title(gamma=-6                  gamma=-5                  gamma=-4)) ytitle("      2-stage Poisson                   2-stage Normal" "     aratio=2        aratio=1        aratio=2        aratio=1") subtitle("")
scatter e b1SCox, mcol(blue) || ///
	line emean zero b1SCox, lcol(blue black) ///
	by(method gamma aratio, legend(off) col(6) note("") ///
		t2title("gamma=-6                                      gamma=-5                                      gamma=-4" "aratio=1                aratio=2                aratio=1                aratio=2                aratio=1                aratio=2")) ///
	ytitle("    log HR, 2-stage Poisson      log HR, 2-stage Normal") ///
	xsize(9) ysize(5) subtitle("") name(discrep_ref,replace) ///
	xlabel(,labsize(large)) ///
	ylabel(,labsize(large)) ///
	xtitle("log HR, 1-stage Cox")

replace method = "2-stage Normal" if method=="2SN"
replace method = "2-stage Poisson" if method=="2SPU"
label var sdb "Standard deviation of point estimates"
scatter e sdb if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sd,replace)
scatter e sdbw if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sdw,replace)

* compare discriminatory ability of sdb and sdbw
gen discrep = abs(e)>=0.05 if !mi(e)
roccomp discrep sdb sdbw, graph name(roccomp,replace) summary plot1opts(ms(i)) plot2(ms(i))


// ANALYSE ESTIMATES OF MU
use twoSPcombine, clear
drop error* tau?*
* drop all 2SP weighted
drop b2SPW* s2SPW*
* drop all 2SP ML uncentred
drop b2SPU_REu s2SPU_REu
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) ///
	method(1SCox 2SN 2SPU 1SWei_REc 2SN_RE 2SPU_REc) true(beta)
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
	lcol(black blue red black blue red) lpat(solid = = dash = =) ///
	name(het,replace)

* summarise MCSEs
table gamma studies if _perf=="pctbias" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Pctbias: Monte Carlo error varies from 1.4 to 2.9 with 5 studies and from 0.7 to 1.3 with 20 studies
table gamma studies if _perf=="empse" & tau>0, stat(min s) stat(max s) nformat(%6.3f)
* Empse: Monte Carlo error < 0.006
table gamma studies if _perf=="cover" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Coverage: Monte Carlo error = 0.4 to 1.0%
table gamma studies if _perf=="power" & tau>0, stat(min s) stat(max s) nformat(%6.1f)
* Power: Monte Carlo error <= 1.1%


// ANALYSE ESTIMATES OF TAU
use twoSPcombine, clear
drop error* b1* b2* s1* s2*
rename tau truetau
* drop all 2SP ML weighted or uncentred
drop tau2SPW* tau2SPU_REu
siman setup, dgm(aratio gamma studies truetau) rep(rep) est(tau) ///
	method(2SN_RE 2SPU_REc) true(truetau)
siman ana, max(100)
siman tab mean if truetau>0, row(truetau gamma aratio studies)
siman nes mean if truetau>0, ///
	legend(pos(3) col(1)) stagger(.05) ///
	$PPT dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(blue red) lpat(dash dash) name(tauhet,replace) yla(0 .1 .2 .3)
table gamma studies if rep<0 & _perf=="mean" & truetau>0, stat(max _se) stat(min _se) nformat(%6.3f)
* Mean: Monte Carlo error < 0.007

/* before using simrecreate,need:
recast double px bx beta
foreach var of varlist px bx beta {
	replace `var' = round(`var',1E-7)
}
*/
