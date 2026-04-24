readme.txt for TwoStagePoisson
==============================

Folder ado
Description
	adofiles for other programs
Contents (packages)
	pllf - compute profile log likelihood
	meta2p (unpublished, ado file only) - two-stage Poisson meta-analysis
	simrun (unpublished) - utility to run simulation studies
	simsum - compute performance measures for simulation studies
	siman - tabulate and graph performance measures for simulation studies

Folder toydata
Description
	data and analyses for the hypothetical data of two studies with grouped Poisson outcome
Inputs
	toydata.dta - the data in Stata format
	toydata.do - analyse by 1-stage and 2-stage Normal. Draw PLLFs with their Normal approximations.
Outputs	
	toydata.log
	toydata_pllf.dta
	toydata_pllf.gph

Folder brcancer
Description
	Analyses of German breast cancer data to illustrate bias due to 2-stage Normal approximation and how 2-stage Poisson analysis can fix it
Inputs
	siloes.do - repeatedly split the data into siloes and analyse by 1-stage and 2-stage Normal. Extract one typical data set and also analyse it by 2-stage Poisson.
Outputs
	siloes.log
	siloes25_11.dta - the one typical data set

Folder simulation
Description
	Simulation program and results
Inputs
	twoSPsimprog.do - defines programs twoSPsimgen to generate data and twoSPsimana to analyse them
	twoSPsimrun.do - uses simrun package to run simulation and compute graphs & tables of performance measures
	twoSPsimrun2.do - as twoSPsimrun.do, but different methods
	twoSPcombine.do - takes both sets of results and combines into twoSPcombine.dta
	twoSPsimresults.do - produce all graphs and tables from twoSPcombine.dta
Outputs
	simrun_results (folder) - estimates datasets and random number generator states
	twoSPsimrun.log - main results
	twoSPcombine.dta - combined data set
	twoSPsimresults.log
