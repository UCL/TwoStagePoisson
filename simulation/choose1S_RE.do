/*
Choose a 1-stage RE analysis
choose1S_RE.do

This program recreates one data set and shows that 
	in CE models, 1S weibull is very close to 1S Cox
	in RE models, mestreg+weibull works nicely
	merlin can in principle fit RE+FPM, but in practice 
		it fails to converge even with CE+exponential 

IW 25feb2026
*/

do twoSPsimrun

simcombine

recast int studies nobs* rep
recast float px aratio gamma bx beta tau cens, force
format px aratio gamma bx beta tau cens %9.0g
de

// let's look at these reps
l rep b1SC s1SC b2SN_REML s2SN_REML tau2SN_REML if studies==20 & tau==float(.3) & gamma==-4 & aratio==1 & rep<=5
/*
       +--------------------------------------------------------------+
       | rep        b1SC       s1SC   b2SN_REML   s2SN_R~L   tau2~EML |
       |--------------------------------------------------------------|
34001. |   1   -.0968821   .0597653   -.0735372   .1134749   .4161574 |
34002. |   2   -.2524757   .0589647   -.2855357   .1072678   .3828048 |
34003. |   3   -.2703566   .0575069   -.2580758   .0837701   .2584515 |
34004. |   4   -.2164178   .0569065   -.2247116   .0713959   .1797116 |
34005. |   5   -.3217182   .0614853    -.389814   .1155417   .4233243 |
       +--------------------------------------------------------------+
*/

// recreate data
simrecreate if studies==20 & tau==float(.3) & gamma==-4 & aratio==1 & rep==1

// check results match
* 1SC
stcox z x, strat(study) nohr

* 2SN_REML
qui metan b se, nograph model(reml)
di r(eff), r(se_eff), sqrt(r(tausq))

// is it ok to replace cox with weibull?
streg z x, strat(study) dist(weibull) nohr
di _b[z], _se[z]

// 1S + REML
* exponential
mestreg z x i.study || study: z, nocons dist(exp) nohr
di _b[z], _se[z], sqrt(_b[/var(z[study])])
* weibull
mestreg z x i.study || study: z, nocons dist(weib) nohr
di _b[z], _se[z], sqrt(_b[/var(z[study])])

// merlin: why it's not OK
* CE, ignore study
stmerlin z x, dist(rp) df(3)

* CE, correctly adjusting for study: none of these converge
forvalues i=1/20 {
	gen Istudy`i'=study==`i'
	local Istudyvars `Istudyvars' Istudy`i'
}
* stmerlin doesn't seem to allow REs, so stick with merlin
* stmerlin z x `Istudyvars' z#M1[study]@1, dist(rp) df(3)

* CE weibull
merlin (t z x `Istudyvars', family(weibull, failure(d)))
* CE FPM
merlin (t z x `Istudyvars', family(rp, df(3) failure(d)))
* RE FPM
merlin (t z x `Istudyvars' z#M1[study]@1, family(rp, df(3) failure(d)))
