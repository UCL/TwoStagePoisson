/* 
simcombine: read and combine all simulation results
*! v0.4 IW 24feb2026
	nosort option keeps sort order of settings frame
		default sorts by inputs
*! v0.3 IW 1oct2025
	new simsetup program that stores settings as char
v0.2 IW 17sep2025
v0.1.1 IW 22jan2021
11jun2020: added genid option
*/

prog def simcombine
version 16
syntax [if] [in], [settings(string) nosort]
if mi("`settings'") local settings simrun_settings

frame `settings' {
	local allthings: char _dta[simrun_allthings] 
	foreach thing of local allthings {
		local `thing' : char _dta[simrun_`thing'] 
	}
}

foreach frame in results data {
	local `frame' simrun_`frame'
}

// END OF PARSING

frame change `settings'
marksample touse
summ `touse', meanonly
local n = r(N)
di as txt "Reading `n' files" _c
local filesused 0
frame `results': clear
forvalues i=1/`n' {
	if !`touse'[`i'] continue
	local postfile = postfilename[`i']
	if `filesused'==0 cap frame `results': use `folder'`postfile'
	else cap frame `results': append using `folder'`postfile'
	if !_rc {
		local ++filesused
		di "." _c
	}
	else di as error "Failed to read `folder'`postfile'"
}
if `filesused'<`n' di as text " only " `filesused' " files found"
else di
frame change `results'
di as text "Results loaded into current frame `results'"
if "`sort'"!="nosort" sort `inputs', stable
if !mi("`genid'") by `inputs': gen `genid' = _n
end
