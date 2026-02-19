/*
meta2p_examples.do
Run the examples in the help file
*/

clear
input study d1 py1 d0 py0
1 40 1000 10 1000
2 60 1000 60 1000
end
gen est = log(d1/py1) - log(d0/py0)
gen se = sqrt(1/d1+1/d0)
metan est se, nograph model(fixed)

/*
local b2sn = r(eff)
qui mvmeta_make, by(study) clear: poisson d z, exp(py) 
mvmeta y S, fixed
assert abs(_b[yz]-`b2sn')<1E-7
*/

meta2p est se, d(d1 d0) 