/*
*! v0.3 IW 1oct2025
	new simsetup program that stores settings as char
NB can be run in frame results or in frame settings
*/

prog def simrecreate
version 16
syntax [if] [in] , [rep(int 0) settings(string) *]
if mi("`settings'") local settings simrun_settings

frame `settings' {
	local allthings: char _dta[simrun_allthings] 
	foreach thing of local allthings {
		local `thing' : char _dta[simrun_`thing'] 
	}
}

foreach frame in settings results data rngstates {
	local `frame' simrun_`frame'
}

// END OF PARSING

qui frame
local frame = r(currentframe) 
if "`frame'"=="`results'" {
	if `rep'>0 di as error "In results frame: rep() ignored"
}
else if "`frame'"=="`settings'" {
	if `rep'==0 {
		di as text "assuming rep(1)"
		local rep 1
	}
}
else {
	di as error "need to run simrecreate in frame `settings' or `results'"
	exit 498
}

marksample touse
qui count if `touse'
if r(N)>1 di as error "Multiple rows identified: using first row"
tempvar id
gen `id'=_n
summ `id' if `touse', meanonly
local row = r(min)
drop `touse' `id'

if "`frame'" == "`results'" {
	* find which row
	local postfilename `name'
	local args
	foreach input of local inputs {
		local inputvalue = `input'[`row']
		local postfilename `postfilename'_`input'`inputvalue'
		local args `args' `input'(`inputvalue')
	}
	local rep = rep[`row']
	local postfilename `postfilename'.dta		
}
else if "`frame'" == "`settings'" {
	local postfilename = postfilename[`row']
	local args = args[`row']
}

frame `rngstates' {
	use `folder'rngstates_`postfilename', clear
	local rngstate = rngstate1[`rep']+rngstate2[`rep']+rngstate3[`rep']
}

// RECREATE THE DATA
di as text "Changing to frame `data'"
frame change `data'
di as text "Running `genprog'"
set rngstate `rngstate'
`genprog', `args' `options'
di as text "Simulated data are now in memory"
end
