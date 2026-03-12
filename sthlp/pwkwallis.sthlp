{smcl}
{* *! version 1.1.1 11jul2025}{...}
{viewerdialog pwkwallis "dialog pwkwallis"}{...}
{vieweralsosee "[R] kwallis" "mansection R kwallis"}{...}
{viewerjumpto "Syntax" "pwkwallis##syntax"}{...}
{viewerjumpto "Description" "pwkwallis##description"}{...}
{viewerjumpto "Examples" "pwkwallis##examples"}{...}
{viewerjumpto "Stored results" "pwkwallis##results"}{...}
{viewerjumpto "Version" "pwkwallis##version"}{...}
{viewerjumpto "Authors" "pwkwallis##authors"}{...}
{viewerjumpto "References" "pwkwallis##references"}{...}
{title:Title}

{phang}
{bf:pwkwallis} {hline 2} Kruskal-Wallis equality-of-populations rank test and non-parametric pairwise 
comparisons across the levels of factor variables


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:pwkwallis} {varname} {ifin}{cmd:,} {opth "by(varlist:groupvar)"} [{cmd:bonferroni}]


{marker description}{...}
{title:Description}

{p 4 4 2}
This command performs Kruskal-Wallis equality-of-populations rank test and non-parametric pairwise 
comparisons across the levels of factor variables. The tests can be adjusted for multiple comparisons 
by Bonferroni's methods using the {cmd:bonferroni} option.

{p 4 4}
In the syntax diagram above, {varname} refers to the variable recording the outcome, 
and {it:{help varlist:groupvar}} refers to the variable denoting the population. 
{opt by()} is required.

{p 4 4}
You can click {dialog pwkwallis:here} to pop up a {dialog pwkwallis:dialog} or type {inp: db pwkwallis}.

{p 4 4}
Execute {cmd: net install pwkwallis, from("https://raw.githubusercontent.com/metodo-leam/stata/master")} for install. 

{p 4 4}
This command uses the {help kwallis} and {help _mtest} Stata commands.


{marker examples}{...}
{title:Examples}

{p 4 4}{stata "use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/AlcoholTR.dta":. use https://raw.githubusercontent.com/metodo-leam/stata/master/dta/AlcoholTR.dta}{p_end}
{p 4 4}{cmd:. pwkwallis TR, by(Alcohol)}{p_end}
{p 4 4}{cmd:. pwkwallis TR, by(Alcohol) bonferroni}{p_end}


{marker results}{...}
{title:Stored results}

{p 4 4}
{cmd:pwkwallis} stores in {cmd:r(results)} a matrix with the contrast, standard error, t and p values (adjusted if {cmd:bonferroni} option
is present) for each pairwise comparison.


{marker version}{...}
{title:Version}

{p 4}
Version 1.1.1 {hline 2} 11 July 2025


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
Dom{c e'}nech JM. Kruskal-Wallis equality-of-populations rank test and non-parametric pairwise comparisons across{break}
the levels of factor variables: User-written command pwkwallis for Stata  [computer program].{break}
V1.1.1. Bellaterra: Universitat Aut{c o'g}noma de Barcelona; 2025.{break}
Available from {browse "https://github.com/metodo-leam/stata"}{p_end}


{marker references}{...}
{title:References}

{p 0 2}
Conover WJ. Practical Nonparametric Statistics. 3th ed. New York: John Wiley and Sons; 1999.{p_end}

{p 0 2}
Dom{c e'}nech JM. Fundamentos de Dise{c n~}o y Estad{c i'}stica. UD 12. Comparaci{c o'}n de varias medias:
An{c a'}lisis de la variancia. 15{c 170} ed. Barcelona: Signo; 2014.{p_end}
