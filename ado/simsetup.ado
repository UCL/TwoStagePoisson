/* 
*! v0.4 IW 26feb2026
	name defaults to current frame name
v0.3 IW 1oct2025
	new simsetup program that stores settings as char
simrun suite
simsetup IW 1oct2025

Frame notes:
	char are frame-specific
	macros are session-specific 
*/

prog def simsetup
version 16
syntax, PROGram(namelist max=2)	/// name of program that generates and analyses data
	[name(name) 	/// file name stub for results and rngstate files
	FOLDer(string)	/// name of folder to contain results and rngstate files (default: simrun_results)
	]

* LEARN NAMES OF DGM PARMS AND OUTPUTS
local nprogs : word count `program'
local genprog : word 1 of `program'
if `nprogs'==2 local anaprog : word 2 of `program'
else local anaprog `genprog'
`genprog' ?
local inputs = r(inputs)
`anaprog' ?
local outputs = r(outputs)

* SORT OUT FRAME NAMES
foreach frame in results data rngstates {
	cap frame create simrun_`frame'
	local `frame' simrun_`frame'
}
qui frame
local settings = r(currentframe)
if mi("`name'") local name `settings'

* SORT OUT FOLDER
if missing("`folder'") local folder simrun_results
* check it exists
local pwd = c(pwd)
cap cd `folder'
if _rc {
	mkdir `folder'
	di as text "Created folder `folder'"
}
qui cd "`pwd'"
if !inlist(substr("`folder'",-1,1),"/","\") local folder `folder'/

// END OF PARSING


// STORE SETTINGS AS CHAR IN FRAME SETTINGS
frame `settings' {
	local allthings name nprogs genprog anaprog inputs outputs folder
	foreach thing of local allthings {
		char _dta[simrun_`thing'] ``thing''
	}
	char _dta[simrun_allthings] `allthings'
}

end	

/* if non-default simsettings() is used then it needs to be passed to all subsequent commands
can't avoid this as simsettings frame is where all locations are stored