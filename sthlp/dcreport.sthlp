{smcl}
{* *! version 1.1.0 12mar2026}{...}
{viewerdialog dcreport "dialog dcreport"}{...}
{vieweralsosee "dc" "help dc"}{...}
{viewerjumpto "Syntax" "dcreport##syntax"}{...}
{viewerjumpto "Description" "dcreport##description"}{...}
{viewerjumpto "Examples" "dcreport##examples"}{...}
{viewerjumpto "Version" "dcreport##version"}{...}
{viewerjumpto "Authors" "dcreport##authors"}{...}
{viewerjumpto "References" "dcreport##references"}{...}
{title:Title}

{phang}
{bf:dcreport} {hline 2} Data Check - Incidence Report


{marker syntax}{...}
{title:Syntax}

{p 4 12 2}
{cmd:dcreport} {varlist} {cmd:,} id(varlist) [idnum]


{marker options}{...}
{synoptset 16 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth id(varlist)}}Identifier variables{p_end}
{synopt :{opt idnum}}Use {bf:_ObsNum} observation number variable to identify cases in the error correction commands{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{p 4 4 2}
Data Check - Incidence Report. 

{p 4 4}
This command uses the auxiliary variables named {bf:_Er_varname} created by the Data Check {cmd:dc} command (see {help dc})
to print an Incidence Report with the problems found. If no variables are specified, all the {bf:_Er_varname}
variables in the dataset are used.

{p 4 4}
An observation number variable named {bf:_ObsNum} must contain the observation number of the original dataset. 

{p 4 4}
You can click {dialog dcreport:here} to pop up a {dialog dcreport:dialog} or type {inp: db dcreport}.

{p 4 4}
Execute {cmd: net install dcreport, from("https://raw.githubusercontent.com/metodo-leam/stata/master")} for install. 


{marker examples}{...}
{title:Examples}
{it:Set data}
{p 4 4}{stata "use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/Salud0.dta":. use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/Salud0.dta}{p_end}
{p 4 4}{cmd:. generate _ObsNum = _n}{p_end}
{it:Check data with dc}
{p 4 4}{cmd:. dc Talla, vl(150/200) id(Id)}{p_end}
{p 4 4}{cmd:. dc Peso, vl(45/120) nd(1) id(Id)}{p_end}
{p 4 4}{cmd:. dc H1 H2 H3, vl(0 1 2) id(Id)}{p_end}

{p 4 4}{cmd:. dcreport, id(Id)}{p_end}
{p 4 4}{cmd:. dcreport _Er_Talla _Er_Peso _Er_H1 _Er_H2 _Er_H3, id(Id) idnum}{p_end}
{p 4 4}{cmd:. dcreport _Er_Talla - _Er_H3, id(Id)}{p_end}


{marker version}{...}
{title:Version}

{p 4}
Version 1.1.0 {hline 2} 12 March 2026


{marker authors}{...}
{title:Authors}

{p 4 4 2}
JM.Dom{c e'}nech{break}
Programmer: R.Sesma{break}
Laboratori d'Estad{c i'}stica Aplicada{break}
Universitat Aut{c o'g}noma de Barcelona{break}
metodo.campus@gmail.com{p_end}


{title:Vancouver reference}

{p 4 6 2}
Dom{c e'}nech JM. Incidence Report: User-written command dcreport for Stata [computer program].{break}
V1.1.0. Bellaterra: Universitat Aut{c o'g}noma de Barcelona; 2026.{break}
Available from {browse "https://github.com/metodo-leam/stata"}{p_end}


{marker references}{...}
{title:References}

{p 0 2}
Bonillo-Mart{c i'}n A. Sistematizaci{c o'}n del proceso de depuraci{c o'}n de los datos en estudios con
seguimientos [Tesis doctoral]. Barcelona: Universitat Aut{c o'g}noma de Barcelona; 2003.{p_end}
