{smcl}
{* *! version 1.0.3 15dec2022}{...}
{viewerdialog intcon "dialog intcon"}{...}
{viewerjumpto "Syntax" "intcon##syntax"}{...}
{viewerjumpto "Description" "intcon##description"}{...}
{viewerjumpto "Examples" "intcon##examples"}{...}
{viewerjumpto "Stored results" "intcon##results"}{...}
{viewerjumpto "Version" "intcon##version"}{...}
{viewerjumpto "Authors" "intcon##authors"}{...}
{viewerjumpto "References" "intcon##references"}{...}
{title:Title}

{phang}
{bf:intcon} {hline 2} Internal Consistency


{marker syntax}{...}
{title:Syntax}

{p 4 12 2}
{cmd:intcon} {varlist} {ifin} {weight}{cmd:, {cont | ord}} [nst]


{marker options}{...}
{synoptset 16 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt cont}}Continuous indexes{p_end}
{synopt :{opt ord}}Ordinal indexes{p_end}
{synopt :{opt nst(string)}}Name of the study (label){p_end}
{synoptline}
{p 4 4 2}{bf:aweight}s, {bf:fweight}s and {bf:pweight}s are allowed. {bf:pweight}s may not be used with continuous indexes.{p_end}


{marker description}{...}
{title:Description}

{p 4 4}
This command computes statistical coefficients to evaluate internal consistency. Cronbach's alpha, Armor's theta and McDonald's omega(t)
are computed for continuous indexes, as well as its equivalents for ordinal indexes. At least two numerical variables must be
specified. At least one {cmd:cont} or {cmd:ord} option is needed.

{p 4 4}
You can click {dialog intcon:here} to pop up a {dialog intcon:dialog} or type {inp: db intcon}.

{p 4 4}
Execute {cmd: net install intcon, from("https://raw.githubusercontent.com/metodo-leam/stata/master")} for install.

{p 4 4}
This command uses the {help alpha}, {help factor}, {help factormat} and {help sem} Stata commands. It also uses the
{stata "net describe alphawgt, from(http://fmwww.bc.edu/RePEc/bocode/a)":alphawgt} and
{stata "net describe polychoric, from(http://staskolenikov.net/stata)":polychoric}
user defined programs.


{marker examples}{...}
{title:Examples}

{p 4 4}{stata "use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/intcon.dta":. use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/intcon.dta}{p_end}
{p 4 4}{cmd:. intcon Item1 Item2 Item3 Item4 Item5 Item6, cont ord}{p_end}

{marker results}{...}
{title:Stored results}

{p 4 4 2}
The command stores the following in {cmd:r()}:

{synoptset 18 tabbed}{...}
{p 4 4}Continuous measures{p_end}
{synopt:{cmd:r(alpha)}}Cronbach's alpha{p_end}
{synopt:{cmd:r(theta)}}Armor's theta{p_end}
{synopt:{cmd:r(omega)}}McDonald's omega(t){p_end}

{p 4 4}Ordinal measures{p_end}
{synopt:{cmd:r(alpha_ord)}}Ordinal alpha{p_end}
{synopt:{cmd:r(theta_ord)}}Ordinal theta{p_end}
{synopt:{cmd:r(omega_ord)}}Ordinal omega{p_end}


{marker version}{...}
{title:Version}

{p 4}
Version 1.0.3 {hline 2} 15 December 2022


{marker authors}{...}
{title:Authors}

{p 4 4 2}
JM.Dom{c e'}nech & JB.Navarro{break}
Programmer: R.Sesma{break}
Laboratori d'Estad{c i'}stica Aplicada{break}
Universitat Aut{c o'g}noma de Barcelona{break}
metodo.campus@gmail.com{p_end}


{title:Vancouver reference}

{p 4 6 2}
Dom{c e'}nech JM, Navarro JB. Internal Consistency: User-written command intcon for Stata [computer program].{break}
V1.0.3. Bellaterra: Universitat Aut{c o'g}noma de Barcelona; 2022.{break}
Available executing from Stata: {break}
net install intcon, from("https://raw.githubusercontent.com/metodo-leam/stata/master"){p_end}


{marker references}{...}
{title:References}

{p 0 2}
Viladrich C, Angulo-Brunet A, Doval E. A journey around alpha and omega to estimate internal consistency reliability. Anal. Psicol. [online]. 2017, 33(3): 755-82. Available from: https://revistas.um.es/analesps/article/view/analesps.33.3.268401
