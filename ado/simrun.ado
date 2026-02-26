/* 
simrun: main program
	List the simulations requested, and run all that need running
*! v0.4 IW 24feb2026
	input variable formats are copied from settings frame to results frame,
		and used in recreating file name
v0.3 IW 1oct2025
	new simsetup program that stores settings as char
v0.2.1 IW 25sep2025
	_dots not dotter
v0.2 IW 17sep2025
	simplified into one program using frames
v0.1.1 IW 22jan2021
	changed variable name _donereps to donereps
	added some documentation
Structure of package:
	User-facing programs:
		simsetup.ado		- set things up
		simrun.ado			- do simulations
		simcombine.ado		- load all simulation results
		simrecreate			- recreate one data set
	Examples:
		simrun.do - define DG parameters and run simple simulation using mysim.ado
		mysim.ado - define DG model, estimands and analyses
Workflow:
	simsetup - load the programs (not the sim settings)
		save settings as chars in frame settings
	simrun - read the simulation settings and perform simulation 
		(i.e. simulate data and analyse data) if needed
	simcombine - combine results
	simrecreate - recreate a data set
Notes 
	donereps=. means there is no file, donereps=0 means an empty file
To do
	If seed is stored in settings file, should it be included in results file name?
	Enable id variable in settings file
	Consider creating a single results file 
		- would avoid float/double difficulties with non-integer dgmvars
*/

prog def simrun
version 16
syntax, [append replace nodots settings(string)]
if mi("`settings'") local settings simrun_settings

frame `settings' {
	local allthings: char _dta[simrun_allthings] 
	foreach thing of local allthings {
		local `thing' : char _dta[simrun_`thing'] 
	}
}

if !mi("`append'") & !mi("`replace'") {
	di as error "Can't use both append and replace"
	exit 198
}

foreach frame in results data rngstates {
	local `frame' simrun_`frame'
}
// END OF PARSING

// LOAD SETTINGS
frame change `settings'
cap confirm var `inputs' reps
if _rc {
	di as error "Error in frame `settings': " _c
	confirm var `inputs' reps
}
cap isid `inputs'
if _rc {
	di as error "{p 0 2}Error in frame `settings': variables `inputs' do not uniquely identify the observations{p_end}"
	exit _rc
}

// CREATE POSTFILENAME, ARGS IF MISSING
cap confirm variable postfilename
if _rc { // postfilename doesn't exist
	confirm new variable args // if postfilename doesn't exist then nor should args
	qui gen postfilename = "`name'"
	qui gen args = ""
	foreach input of local inputs {
		local format : format `input'
		qui replace postfilename = postfilename + "_" + "`input'" + string(`input',"`format'")
		qui replace args = args + " `input'(" + string(`input',"`format'") + ")"
	}
	qui replace postfilename = postfilename + ".dta"
}

// RECREATE DONEREPS
cap drop donereps
qui gen donereps = .
forvalues i=1/`=_N' {
	local postfilename = postfilename[`i']
	cap confirm file `folder'`postfilename'
	if !_rc { // file exists
		if !mi("`append'") { // append to existing files
			qui desc using `folder'`postfilename'
			qui replace donereps = r(N) in `i'
		}
		else if !mi("`replace'") { // replace existing files
			erase `folder'`postfilename'
			di as text "Deleted file: " as result "`postfilename'"
		}
		else {
			di as error "`postfilename' already exists: use replace or append option"
			exit 498
		}
	}
	else { // create files
		confirm new file `folder'`postfilename'
	}
}

// DISPLAY CURRENT STATE
di as text "Current state of simulation:"
l `inputs' reps donereps 

// FIND SOMETHING TO DO
cap assert max(0,donereps)>=reps 
if !_rc {
	di as text "Nothing to do"
	exit
}

forvalues i=1/`=_N' {
	frame change `settings'
	cap frame drop `results'
	cap frame drop `rngstates'
	local prevreps = max(0,donereps[`i']) 
	local reps = max(0, reps[`i'])
	if `prevreps' >= `reps' continue // this line is done
	local postfilename = postfilename[`i']
	local args = args[`i']
	local todo = `reps' - `prevreps'
	local prevrepsplus1 = `prevreps'+1
	di as text _new "Running repetitions " as result `prevrepsplus1' as text " to " as result `reps' as text " of row " as result `i' as text ":"
	if `prevreps'==0 {
		* frame create `results' double(`inputs') int(rep) `outputs' 
		* above loses input variable formats
		frame copy `settings' `results'
		frame `results' {
			qui drop in 1/l
			keep `inputs'
			order `inputs'
			gen int rep = .
			foreach var in `outputs' {
				gen float `var' = .
			}
		}
		* non-integer inputs as double ensures values match in simrecreate
		frame create `rngstates' int(rep) str2000(rngstate1 rngstate2 rngstate3) 
		cap confirm var seed
		if !_rc {
			local seed = seed[`i']
			set seed `seed'
			di as text "Setting seed to " as result `seed'
		}
		else di as error "seed not set: results are dependent on #reps"
		local create create
	}
	else { // ALREADY HAVE SOME RESULTS
		frame create `results'
		frame create `rngstates' 
		frame `results': use `folder'`postfilename', clear
		frame `rngstates': use `folder'rngstates_`postfilename', clear
		frame `results': local nsimrun = _N
		frame `rngstates': local nrngstates = _N
		if `nrngstates'!=`nsimrun'+1 di as error "{p 0 2}`results' and `rngstates' file have non-matching #obs `nsimrun' `nrngstates'{p_end}"
		frame `rngstates': local rngstate = rngstate1[_N]+rngstate2[_N]+rngstate3[_N]
		frame `rngstates': qui drop in l
		set rngstate `rngstate'
		local create append to
	}
	local postin
	foreach parm of local inputs {
		local parmvalue = `parm'[`i']
		local postin `postin' (`parmvalue')
	}

	// RUN THE SIMULATION
	cap frame create simdata
	frame change simdata
	forvalues thisrep=0/`todo' {
		local rep = `prevreps'+`thisrep'+1
		* post RNG states
		local rngstate1 = substr(c(rngstate),1,2000)
		local rngstate2 = substr(c(rngstate),2001,2000)
		local rngstate3 = substr(c(rngstate),4001,.)
		frame post `rngstates' (`rep') ("`rngstate1'") ("`rngstate2'") ("`rngstate3'")
		if mi("`dots'") _dots `thisrep' result
		if `thisrep'==`todo' continue, break

		* generate and analyse
		cap `genprog', `args' 
		if _rc==1 {
			local userbreak 1
			continue, break // responds to user break
		}
		else if _rc {
			di as error "`genprog' exited with error " _rc 
			exit _rc
		}
		if `nprogs'==2 {
			cap `anaprog'
			if _rc==1 {
				local userbreak 1
				continue, break // responds to user break
			}
			else if _rc {
				di as error "`anaprog' exited with error " _rc 
				exit _rc
			}
		}
		local donereps `rep'
		
		* post results
		local postout (`rep')
		foreach parm of local outputs {
			local parmvalue = r(`parm')
			local postout `postout' (`parmvalue')
		}
		frame post `results' `postin' `postout'
		
	}

	frame `results': qui save `folder'`postfilename', replace
	frame `rngstates': qui save `folder'rngstates_`postfilename', replace
	frame `settings': qui replace donereps = `donereps' in `i'
	if "`userbreak'"=="1" {
		di as error "Break after repetition `donereps'"
		exit 1
	}
	if mi("`dots'") di
	di as text "Simulation complete, file saved"

	* check results
	frame `results' {
		local maxnmiss 0
		local n=_N
		local allmissvars
		local partmissvars
		foreach output of local outputs {
			qui count if mi(`output')
			if r(N)==_N {
				local allmissvars `allmissvars' `output'
			}
			else if r(N)>0 {
				local partmissvars `partmissvars' `output' (`=r(N)')
			}
		}
		if !mi("`allmissvars'") {
			di as error "{p 0 2}Wholly missing outputs: `allmissvars'{p_end}"
		}
		if !mi("`partmissvars'") {
			di as error "{p 0 2}Incomplete outputs (# missing): `partmissvars'{p_end}"
		}
	}

} // END OF LOOP
di as text _n "All simulations completed"
frame change `settings' // return to original frame
end	
