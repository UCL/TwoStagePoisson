/*
twoSPsimresults.do
Simulation study to evaluate two-stage Poisson (2SP) method for meta-analysis of TTE data
Main results for paper
*/

set linesize 120

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
mvpatterns bCox* bWei*
mvpatterns bN* bP*


// COMPARE REFERENCE METHODS
use twoSPcombine, clear
gen absdiff=abs(bWei - bCox)
su absdiff
* Wei and Cox never differ by more than 0.004

drop error* tauCox* tauWei* tauN* tauP* 
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(Cox Cox_i Wei) true(beta)
siman ana, max(15)
siman nes pctbias, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black green purple) lpat(solid = =) ///
	name(compareC,replace) yla(-15(5)5)
siman nes modelse, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black green purple) lpat(solid = =) ///
	name(compareC,replace) yla(0 .1 .2 .3)
* the 3 methods have indistinguishable bias and modelSE
table studies method if rep==-4, stat(min s) stat(max s)
table studies method if rep==-6, stat(min s) stat(max s)


// COMPARE UNWEIGHTED AND WEIGHTED P METHODS
use twoSPcombine, clear
drop error* tauCox* tauWei* tauN* tauP*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(PU PW PU_REc PW_REc) true(beta)
siman ana
siman nes pctbias ///
	if inlist(method,"PU":methodlabel,"PW":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingCE,replace) yla(-15(5)5)
siman nes modelse ///
	if inlist(method,"PU":methodlabel,"PW":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingCE,replace) yla(.1 .2 .3)
siman nes pctbias ///
	if inlist(method,"PU_REc":methodlabel,"PW_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingRE,replace) yla(-15(5)5)
siman nes modelse ///
	if inlist(method,"PU_REc":methodlabel,"PW_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(weightingRE,replace) yla(.1 .2 .3 .4)
table gamma studies if rep<0 & _perf=="pctbias", stat(max s) stat(min s) nformat(%6.3f)
* pctbias: Monte Carlo error varies from 0.9 to 2.7 with 5 studies and from 0.4 to 1.3 with 20 studies
* modelse: Monte Carlo error < 0.003


// COMPARE UNCENTRED AND CENTRED METHODS
use twoSPcombine, clear
drop error* tauCox* tauWei* tauN* tauP*
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) method(Wei_REc Wei_RE PU_REc PU_RE) true(beta)
siman ana
siman nes pctbias ///
	if inlist(method,"PU_RE":methodlabel,"PU_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(centringPU,replace) yla(-15(5)5)
siman nes cover ///
	if inlist(method,"PU_RE":methodlabel,"PU_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(red =) lpat(solid dash) ///
	name(centringPU,replace) yla(80(5)100)
siman nes pctbias ///
	if inlist(method,"Wei_RE":methodlabel,"Wei_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black =) lpat(solid dash) ///
	name(centringWei,replace) yla(-15(5)5)
siman nes cover ///
	if inlist(method,"Wei_RE":methodlabel,"Wei_REc":methodlabel), ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 = = =) dglwidth(*2) ///
	lcol(black =) lpat(solid dash) ///
	name(centringWei,replace) yla(80(5)100)
table studies method if rep==-4, stat(min s) stat(max s)
* pctbias: Monte Carlo error varies from 0.9 to 2.9 with 5 studies and from 0.4 to 1.3 with 20 studies
table studies method if rep==-13, stat(min s) stat(max s)
* cover: Monte Carlo error varies from 0.4 to 0.9%.


// SCATTER PLOT OF DIFFERENCES FROM C: HOM, CE METHODS
use twoSPcombine, clear
drop error* tauCox* tauWei* tauN* tauP* sCox* sWei* sN* sP*
foreach method in N PU {
	gen e`method' = b`method' - bCox
}
reshape long e, i(studies aratio gamma tau rep) j(method) string
label var e "Difference from Cox"
label var bCox "Cox"
egen emean = mean(e), by(aratio method)
gen zero=0
recast double tau
replace tau = round(tau,.1)
keep if studies==5 & tau==0
*scatter e bC, mcol(blue) || line emean zero bC, lcol(blue black) by(method aratio gamma, legend(off) col(3) note("") t2title(gamma=-6                  gamma=-5                  gamma=-4)) ytitle("      2-stage Poisson                   2-stage Normal" "     aratio=2        aratio=1        aratio=2        aratio=1") subtitle("")
scatter e bCox, mcol(blue) || ///
	line emean zero bCox, lcol(blue black) ///
	by(method gamma aratio, legend(off) col(6) note("") ///
		t2title("gamma=-6                                      gamma=-5                                      gamma=-4" "aratio=1                aratio=2                aratio=1                aratio=2                aratio=1                aratio=2")) ///
	ytitle("    log HR, Poisson approx     log HR, Normal approx") ///
	xsize(9) ysize(5) subtitle("") name(discrep_ref,replace) ///
	xlabel(,labsize(large)) ///
	ylabel(,labsize(large)) ///
	xtitle("log HR, Cox")

replace method = "Normal approx" if method=="N"
replace method = "Poisson approx" if method=="PU"
label var sdb "Standard deviation of point estimates"
scatter e sdb if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sd,replace)
scatter e sdbw if gamma==-6 & aratio==1, by(method, legend(off) note("")) $PPT name(discrep_sdw,replace)

* compare discriminatory ability of sdb and sdbw
gen discrep = abs(e)>=0.05 if !mi(e)
roccomp discrep sdb sdbw, graph name(roccomp,replace) summary plot1opts(ms(i)) plot2(ms(i))


// ANALYSE ESTIMATES OF MU
use twoSPcombine, clear
drop error* tauCox* tauWei* tauN* tauP*
* drop all P weighted
drop bPW* sPW*
* drop all P ML uncentred
drop bPU_RE sPU_RE
siman setup, dgm(aratio gamma studies tau) rep(rep) est(b) se(s) ///
	method(Cox N PU Wei_REc N_RE PU_REc) true(beta)
siman ana

* bias & power, homogeneous & CE methods
siman tab pctbias empse cover power if tau==0 & method<=3, col(method gamma studies) row(_perf aratio) nformat(%6.1f)
siman nes pctbias if tau==0 & method<=3, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) yla(-10 0 10 20) ///
	lcol(black blue red =) lpat(solid = = dash) name(homCE,replace)
siman nes empse if tau==0 & method<=3, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) yla(0 .1 .2 .3 .4) ///
	lcol(black blue red =) lpat(solid = = dash) name(homCE,replace)
siman nes cover if tau==0 & method<=3, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) yla(94/98) ///
	lcol(black blue red =) lpat(solid = = dash) name(homCE,replace)
siman nes power if tau==0 & method<=3, ///
	legend(pos(3) col(1)) stagger(.05) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) yla(0 50 100) ///
	lcol(black blue red =) lpat(solid = = dash) name(homCE,replace)
* summarise MCSEs
table gamma studies if rep<0 & _perf=="pctbias" & tau==0 & method<=3, stat(max s) stat(min s) nformat(%6.1f)
* pctbias: MCSE varies from 0.9 to 2.7 with 5 studies and from 0.4 to 1.2 with 20 studies
* Empse: MCSE < 0.006
* Coverage: MCSE = 0.3 to 0.5%
* power: MCSE <= 1.1%

* bias, heterogeneous & RE methods
siman tab pctbias empse cover power if tau>0, col(method gamma studies) row(_perf aratio) nformat(%6.1f)
siman nes pctbias if tau>0, ///
	legend(pos(3) col(1)) stagger(.08) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red black blue red) lpat(solid = = dash = =) ///
	name(het,replace) yla(-30 -20 -10 0) debug
siman nes empse if tau>0, ///
	legend(pos(3) col(1)) stagger(.08) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red black blue red) lpat(solid = = dash = =) ///
	name(het,replace) 
siman nes cover if tau>0, ///
	legend(pos(3) col(1)) stagger(.08) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red black blue red) lpat(solid = = dash = =) ///
	name(het,replace) yla(70 80 90 100)
siman nes pctbias empse cover power if tau>0, ///
	legend(pos(3) col(1)) stagger(.08) $PPT ///
	dglabsize(small) lw(*2 ...) dglwidth(*2) ///
	lcol(black blue red black blue red) lpat(solid = = dash = =) ///
	name(het,replace) yla(0 50 100)

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
drop error* bCox* bWei* bN* bP* sCox* sWei* sN* sP*
rename tau truetau
* drop all P ML weighted or uncentred
drop tauPW* tauPU_RE
siman setup, dgm(aratio gamma studies truetau) rep(rep) est(tau) ///
	method(N_RE PU_REc) true(truetau)
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
