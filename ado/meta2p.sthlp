{smcl}
{* *! version 1.0 17 Feb 2026}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "meta2p##syntax"}{...}
{viewerjumpto "Description" "meta2p##description"}{...}
{viewerjumpto "Options" "meta2p##options"}{...}
{viewerjumpto "Remarks" "meta2p##remarks"}{...}
{viewerjumpto "Examples" "meta2p##examples"}{...}
{title:Title}
{phang}
{bf:meta2p} {hline 2} Two-stage Poisson approximation for Meta-Analysis

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:meta2p}
est se
[{help if}]
[{help in}]
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt d(d1 d2)}}Names of variables containing event counts in groups 1 and 2.

{syntab:Optional}
{synopt:{opt s:tudy(varname)}}Study identifier{p_end}
{synopt:{opt py:ears(p1 p2)}}Names of variables containing person-years in groups 1 and 2. Only used if single-zero studies are found.{p_end}
{synopt:{opt wt}}Use weights to make the Poisson data exactly reproduce the observed standard errors.{p_end}
{synopt:{opt re}}Fit random-effects model. Default is a common-effect model.{p_end}
{synopt:{opt cen:tre}}Centre the treatment covariate.{p_end}
{synopt:{opt verb:ose}}Show the output of the Poisson regression.{p_end}
{synopt:{opt pause}}Pause before the Poisson regression. Useful for understanding errors.{p_end}
{synopt:{opt list}}Summarise Poisson data before fitting Poisson model.{p_end}
{synopt:{opt irr}}Report exponentiated coefficient (incidence rate ratio). Same as eform.{p_end}
{synopt:{opt eform}}Report exponentiated coefficient (incidence rate ratio). Same as irr.{p_end}
{synopt:{opt poisson:options(string)}}Options for the Poisson regression. Rarely used.{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:meta2p} performs meta-analysis of estimated treatment effects in time-to-event data, 
avoiding the Normal approximation. Instead, a Poisson likelihood is approximated by creating pseudo-data with the observed numbers 
of events but artificial numbers of person-years such that the point estimate is reproduced exactly. A Poisson regression is then 
fitted to the pseudo-data. Optionally, weights are applied to reproduce the standard errors exactly.

{pstd} 
{it:est} must be the log hazard ratio (or log rate ratio) comparing group 1 with group 2, and {it:se} must be its standard error.

{marker options}{...}
{title:Options}
{dlgtab:Main}

{phang}{opt d(d1 d2)}  Names of variables containing event counts in groups 1 and 2. Required.

{phang}{opt py:ears(p1 p2)}  Names of variables containing person-years in groups 1 and 2. Only used if single-zero studies are found.{p_end}

{phang}{opt wt}  Use weights to make the Poisson data exactly reproduce the observed standard errors.{p_end}

{phang}{opt re}  Fit random-effects model. Default is a common-effect model.{p_end}

{phang}{opt cen:tre}  Centre the treatment covariate.  This gives better results for the between-study heterogeneity variance, when fitting a random-effects model.{p_end}

{phang}{opt irr}  Report exponentiated coefficient (incidence rate ratio). Same as eform.{p_end}

{dlgtab:Other}

{phang}{opt s:tudy(varname)}  Study identifier. Has no effect.{p_end}

{phang}{opt verb:ose}  Show the output of the Poisson regression.{p_end}

{phang}{opt pause}  Pause before the Poisson regression. Useful for understanding errors.{p_end}

{phang}{opt list}  Summarise Poisson data before fitting Poisson model.{p_end}

{phang}{opt eform}  Report exponentiated coefficient (incidence rate ratio). Same as irr.{p_end}

{phang}{opt poisson:options(string)}  Options for the Poisson regression. Rarely used.{p_end}


{marker examples}{...}
{title:Examples}

{phang}Input two-group data with event outcome

{phang}. {stata "clear"}{p_end}
{phang}. {stata "input study d1 py1 d0 py0"}{p_end}
{phang}{space 4}{stata "1 40 1000 10 1000"}{p_end}
{phang}{space 4}{stata "2 60 1000 60 1000"}{p_end}
{phang}{space 4}{stata "end"}{p_end}

{phang}Two-stage common-effect Normal-approximation meta-analysis

{phang}. {stata "gen est = log(d1/py1) - log(d0/py0)"}{p_end}
{phang}. {stata "gen se = sqrt(1/d1+1/d0)"}{p_end}
{phang}. {stata "metan est se, nograph model(fixed)"}{p_end}

{phang}Two-stage common-effect Poisson-approximation meta-analysis

{phang}. {stata "meta2p est se, d(d1 d0)"}{p_end}

{phang}Weighted version (gives same answers here)

{phang}. {stata "meta2p est se, d(d1 d0) wt"}{p_end}


{title:Stored results}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(eff)}}pooled estimate {p_end}
{synopt:{cmd:r(se_eff)}}standard error of pooled estimate  {p_end}
{synopt:{cmd:r(tausq)}}estimate of tau-squared (if re option used)  {p_end}
{synopt:{cmd:r(error)}}flag of any problem in the Poisson analysis  {p_end}


{title:Author}
{pstd}Ian White, MRC Clinical Trials Unit at UCL, London, UK.{break}
Email: {browse "mailto:ian.white@ucl.ac.uk":Ian White}



