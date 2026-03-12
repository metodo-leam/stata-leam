*! version 1.2.4  15sep2025 JM. Domenech, R. Sesma
/*
Association Measures: FREQUENCY DATA
*/

program define sta__freq, rclass
	version 12
	syntax [anything], [Data(string) ST(string) Level(numlist max=1 >50 <100) Wilson  	/*
	*/	 Exact WAld Pearson MH R(numlist max=1 >0 <1) PE(numlist max=1 >0 <1) DETail 	/*
	*/	 notables NNT(numlist integer max=1 >=0 <=1) Zero(string) rare nst(string)		/*
	*/   _incall _data(name) _valid(integer 0) _res(varname) _exp(varname) _touse(varname)]

	* default level: 95
	if ("`level'"=="") local level 95
	* check st; cs by default
	if ("`st'"=="") local st = "cs"
	if ("`st'"!="cs" & "`st'"!="co" & "`st'"!="ex" & "`st'"!="cc") print_error "st() invalid -- invalid value"
	* check confidence interval method; wilson by default
	local ci_type = trim("`wilson' `exact' `wald'")
	if ("`ci_type'"=="") local ci_type = "wilson"
	if (`: word count("`ci_type'")'>1) print_error "only one of wilson, exact, wald options is allowed"
	local ci_lb = upper(substr("`ci_type'",1,1)) + substr("`ci_type'",2,.)
	* check chi2 method; pearson by default
	local chi2_type = trim("`pearson' `mh'")
	if (`: word count("`chi2_type'")'>1) print_error "pearson and mh options are incompatible -- choose one"
	if ("`chi2_type'"=="") local chi2_type = "pearson"		// pearson by default
	local chi2_type = cond("`chi2_type'"=="pearson",1,2)
	* check pe, r, rare, nnt
	if ("`pe'"=="") local pe = 0
	if ("`r'"=="") local r = 0
	if ("`st'"=="cc" & `r'!=0 & `pe'!=0) print_error "options r and pe are incompatible -- choose one"
	if ("`st'"!="cc" & `r'!=0) print_error "option r is only available for case-control studies"
	if ("`st'"!="cc" & "`rare'"!="") print_error "option rare is only available for case-control studies"
	if ("`st'"=="cc" & "`rare'"=="" & `r'==0 & `pe'==0 & "`detail'"!="") print_error "for not rare disease, detail results require r or pe"
	if ("`st'"!="ex" & "`st'"!="co" & "`nnt'"!="") print_error "option nnt is only available for cohort and experimental studies"
	if ("`nnt'"=="") local nnt = -1
	local det = ("`detail'"!="")
	local rare_disease = ("`rare'"!="")
	* check zero correction: none by default
	if ("`zero'"=="") local zero = "n"
	if ("`zero'"!="n" & "`zero'"!="c" & "`zero'"!="p" & "`zero'"!="r") print_error "zero() invalid -- invalid value"
	
	if ("`relatsymm'"!="") print_error "relatsymm option only makes sense with paired data"

	*** GET DATA
	* a0 a1 b0 b1
	if ("`_incall'"!="") {
		local a0 = `_data'[2,1]
		local a1 = `_data'[2,2]
		local b0 = `_data'[1,1]
		local b1 = `_data'[1,2]
	}
	else {
		tokenize `anything'
		confirm number `1'
		confirm number `2'
		confirm number `3'
		confirm number `4'
		local a1 = `1'
		local a0 = `2'
		local b1 = `3'
		local b0 = `4'
	}
	* m1 m0 n0 n1 n
	local m1 = `a0' + `a1'
	local m0 = `b0' + `b1'
	local n0 = `a0' + `b0'
	local n1 = `a1' + `b1'
	local n = `a0' + `b0' + `a1' + `b1'
	* data matrix
	tempname D
	matrix `D' = (`a0',`a1',`m1' \ `b0',`b1',`m0' \ `n0',`n1',`n')
		
	*** GET RESULTS
	sta__utils is_dec, data(`D')			// is there decimal values?
	local ldec = r(dec)
		
	* get chi2 results before any continuity correction
	tempname C
	sta__utils get_chi2_fisher, data(`D') type(`chi2_type') st(`st') dec(`ldec')
	matrix define `C' = r(chi2)
	
	local lzero 0
	local z " "
	if ("`zero'"!="n") {
		sta__utils is_zero, data(`D')		// is there some zero value?
		local lzero = r(zero)
		if (`lzero') {
			local z "*"
			sta__utils zero_correction `zero', d(`D') st(`st')
			local a0 = `D'[1,1]
			local a1 = `D'[1,2]
			local m1 = `D'[1,3]
			local b0 = `D'[2,1]
			local b1 = `D'[2,2]
			local m0 = `D'[2,3]
			local n0 = `D'[3,1]
			local n1 = `D'[3,2]
			local n = `D'[3,3]
		}
	}
	local ldec = (`ldec' | `lzero')
	
	* get proportions
	tempname P
	sta__utils get_proportions, d(`D') level(`level') method(`ci_type')
	matrix `P' = r(p)

	* get estimations
	tempname R
	if ("`st'"!="cc") sta__utils get_cor_results, d(`D') st(`st') level(`level') nnt(`nnt') pe(`pe') dec(`ldec') detail(`det')
	else sta__utils get_cc_results, d(`D') st(`st') level(`level') r(`r') pe(`pe') rare(`rare_disease') dec(`ldec') detail(`det')
	matrix define `R' = r(results)


	*** PRINT RESULTS
	* print title
	if ("`st'"=="cs") di as res "ASSOCIATION MEASURES: CROSS-SECTIONAL STUDY"
	if ("`st'"=="co") di as res "ASSOCIATION MEASURES: COHORT STUDY"
	if ("`st'"=="ex") di as res "ASSOCIATION MEASURES: EXPERIMENTAL STUDY"
	if ("`st'"=="cc") di as res "ASSOCIATION MEASURES: CASE-CONTROL STUDY"
	if ("`nst'"!="") di as txt "{bf:STUDY:} `nst'"
	* print valid/total observations
	if ("`_incall'"!="") {
		qui count if `_touse'			// count number of total cases
		local total = r(N)
		di as txt "Valid observations: " as res `_valid' as txt " (" as res %5.1f 100*`_valid'/`total' as txt "%)"
		di as txt "Total observations: " as res `total'
	}
	
	*** PRINT TABLES
	if ("`tables'"=="") {
		* get column / row labels
		if ("`_incall'"!="") {
			* get variable & value labels for exposure var
			sta__utils get_var_labels `_exp', abb_name(21) abb_lbl(10)
			local col_header = r(vname)
			local col_lb1 = r(lbl0)
			local col_lb2 = r(lbl1)
			* get variable & value labels for response var
			sta__utils get_var_labels `_res', abb_name(17) abb_lbl(17)
			local row_header = r(vname)
			local row_lb1 = r(lbl1)
			local row_lb2 = r(lbl0)
		}
		else {
			* default labels for immediate calls
			local row_header
			local col_header
			local col_lb1 = cond("`st'"=="ex","Treat.0 ","Unexposed ")
			local col_lb2 = cond("`st'"=="ex","Treat.1 ","Exposed ")
			local row_lb1 = cond("`st'"=="cc","Cases","Events")
			local row_lb2 = cond("`st'"=="cc","Controls","NonEvents")
		}
		
		di
		if ("`st'"=="cs" | "`st'"=="cc") {
			* cross-sectional & case control show exposed proportion
			* header
			if ("`_incall'"!="") di as txt _col(19) "{c |}{bf:{center 21:`col_header'}}" _col(41) "{c |}" _c
			di as txt _col(52) "{c |}  Exposed proportion"
			di as txt "{bf:{ralign 17:`row_header'}} {c |}{ralign 10:`col_lb1'}{c |}{ralign 10:`col_lb2'}" /*
			*/ 			"{c |}{ralign 9:TOTAL} {c |}  `level'% CI (`ci_lb')"
			di "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}{c +}{hline 23}"
			* body
			di as txt "{ralign 17:`row_lb1'} {c |} " as res %8.0g `a0' as txt " {c |} " as res %8.0g `a1' /*
			*/	as txt " {c |} " as res %8.0g `m1' as txt " {c |} " as res %8.0g `P'[1,1] "`z'"
			di as txt _col(19) "{c |}" _col(30) "{c |}" _col(41) "{c |}" _col(52) "{c |} " /*
			*/  as txt "(" as res %-8.0g `P'[1,2] as txt " to " as res %8.0g `P'[1,3] as txt ")"
			di as txt "{ralign 17:`row_lb2'} {c |} " as res %8.0g `b0' as txt " {c |} " as res %8.0g `b1' /*
			*/	as txt " {c |} " as res %8.0g `m0' as txt " {c |} " as res %8.0g `P'[2,1] "`z'"
			di as txt _col(19) "{c |}" _col(30) "{c |}" _col(41) "{c |}" _col(52) "{c |} " /*
			*/  as txt "(" as res %-8.0g `P'[2,2] as txt " to " as res %8.0g `P'[2,3] as txt ")"
			di "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}{c +}{hline 23}"
			di as txt "{ralign 17:TOTAL} {c |} " as res %8.0g `n0' as txt " {c |} " as res %8.0g `n1' /*
			*/	as txt " {c |} " as res %8.0g `n' as txt " {c |}" _c
			if ("`st'"=="cs") {
				di as res " " %8.0g `P'[3,1] "`z'"
				di as txt _col(19) "{c |}" _col(30) "{c |}" _col(41) "{c |}" _col(52) "{c |} " /*
				*/ as txt "(" as res %-8.0g `P'[3,2] as txt " to " as res %8.0g `P'[3,3] as txt ")"
			}
			if ("`st'"=="cc" & (`r'!=0 | `pe'!=0)) {
				* special case: case-control + additional information
				sta__utils get_risk_odds, data(`D') r(`r') pe(`pe')
				di as res "  " %-8.0g r(pext) " (`r(text)')"
				local rpext = r(pext)		//save for later: stored results need that value
				di _col(19) "{c |}" _col(30) "{c |}" _col(41) "{c |}" _col(52) "{c |}"
				di as txt "{ralign 17:Odds} {c |} " as res %8.0g r(o0) "`z'{c |} " %8.0g r(o1) "`z'{c |} " %8.0g r(odds) "`z'{c |}"
				di as txt "{ralign 17:Proportions} {c |} " as res %8.0g r(r0) "`z'{c |} " %8.0g r(r1) "`z'{c |} " %8.0g r(re) "`z'{c |}"
				di
			}
			if ("`st'"=="cc" & (`r'==0 & `pe'==0)) di
		}
		if ("`st'"=="co" | "`st'"=="ex") {
			* cohort & experimental don't show exposed proportion
			* header
			if ("`_incall'"!="") di as txt _col(19) "{c |}{bf:{center 21:`col_header'}}" _col(41) "{c |}"
			di as txt "{bf:{ralign 17:`row_header'}} {c |}{ralign 10:`col_lb1'}{c |}{ralign 10:`col_lb2'}{c |}{ralign 9:TOTAL}"
			di "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}"
			* body
			di as txt "{ralign 17:`row_lb1'} {c |} " as res %8.0g `a0' as txt " {c |} " as res %8.0g `a1' as txt " {c |} " as res %8.0g `m1'
			di as txt "{ralign 17:`row_lb2'} {c |} " as res %8.0g `b0' as txt " {c |} " as res %8.0g `b1' as txt " {c |} " as res %8.0g `m0'
			di "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}"
			di as txt "{ralign 17:TOTAL} {c |} " as res %8.0g `n0' as txt " {c |} " as res %8.0g `n1' as txt " {c |} " as res %8.0g `n'
			di as txt _col(19) "{c |}" _col(30) "{c |}" _col(41) "{c |}"
		}
		* footer
		if ("`st'"!="cc") {
			local c = cond("`st'"=="cs","Prevalence","Risk")
			di as txt "{ralign 17:`c'} {c |} " as res %8.0g `P'[4,1] as txt "`z'{c |} " as res %8.0g `P'[5,1] as txt "`z'{c |} " _c
			if ("`st'"=="cs") di as res %8.0g `P'[6,1] as txt "`z'{c |}" _c
			di as txt _n %4.0g `level' "% CI" _col(13) "Lower {c |} " as res %8.0g `P'[4,2] as txt " {c |} " as res %8.0g `P'[5,2] as txt " {c |} " _c
			if ("`st'"=="cs") di as res %8.0g `P'[6,2] as txt " {c |}" _c
			di as txt _n " (`ci_lb')" _col(13) "Upper {c |} " as res %8.0g `P'[4,3] as txt " {c |} " as res %8.0g `P'[5,3] as txt " {c |} " _c
			if ("`st'"=="cs") di as res %8.0g `P'[6,3] as txt " {c |}" _c
			di
		}
	}
	
	*** PRINT ESTIMATIONS
	local nnt_h ""
	if (("`st'"=="ex" | "`st'"=="co") & `nnt'>=0) local nnt_h = cond(`nnt'==1,"Beneficial events","Harmful events")
	local se_h = "Std. Err."
	if ("`st'"=="cc") local se_h = "SE(lnOR)"
	sta__utils get_note			// get ASCII value for note Recommended CI
	local note = r(note)

	di
	di as txt "{hline 22}{c TT}{hline 12}{c TT}{hline 31}"
	di as txt "`nnt_h'{col 23}{c |}{col 26}Estimate{col 36}{c |}{ralign 10:`se_h'}  [`level'% Conf. Interval]"
	di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
	
	if ("`st'"!="cc") {
		** CROSS-SECTIONAL, COHORT, EXPERIMENTAL studies
		* Prevalence / Risk difference PD/RD - R[1,.] (CI Newcombe) + R[2,.] (CI Wald)
		local lbl = cond("`st'"=="cs","Prev. Diff. (PD)","Risk Diff. (RD)")
		di as txt "{ralign 21: `lbl'} {c |} " as res %9.0g `R'[1,1] as txt "`z' {c |} " /*
		*/	_col(49) as res %9.0g `R'[1,2] " " %9.0g `R'[1,3] " Newcombe{c `note'}"
		di as txt _col(23) "{c |} " _col(36) "{c |} " as res %9.0g `R'[2,4] "  " as res %9.0g `R'[2,2] " " %9.0g `R'[2,3] as txt " Wald"
		di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
		
		if (`det'==1) {
			* PRINT DETAIL RESULTS (if asked)
			* Proportion of exposed in the population - R[7,1]
			di as txt "{ralign 21: Prop. exposed pop.} {c |} " as res %9.0g `R'[7,1] "`z' {c |} " cond(`pe'==0,"(estimated)","(external)")
			
			* RDp/PDp: Risk / Prevalence difference in the population - R[8,.]
			local lbl = cond("`st'"=="cs","Prev. Diff. pop.","Risk Diff. pop.")
			di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[8,1] "`z' {c |} " _c
			if (`pe'==0) di as res %9.0g `R'[8,4] "  " _c
			else di _col(49) _c
			di as res %9.0g `R'[8,2] " " %9.0g `R'[8,3]
			
			* AF/PFe: Attributable (R[9,4]=1) / Preventable (R[9,4]=2) fraction in exposed - R[9,.]
			local lbl = cond(`R'[9,4]==1,"Attr. Frac. exp.","Prev. Frac. exp.")
			di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[9,1] "`z' {c |} " _col(49) _c
			if (`R'[9,2]<. | `R'[9,3]<.) di as res  %9.0g `R'[9,2] " " %9.0g `R'[9,3]
			else di
			
			* AF/PFp: Attributable (R[9,4]=1) / Preventable (R[9,4]=2) fraction in population - R[10,.]
			local lbl = cond(`R'[9,4]==1,"Attr. Frac. pop.","Prev. Frac. pop.")
			di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[10,1] "`z' {c |} " _c
			if (`R'[10,4]<.) di as res %9.0g `R'[10,4] "  " _c
			di _col(49) _c
			if (`R'[10,2]<. | `R'[10,3]<.) di as res  %9.0g `R'[10,2] " " %9.0g `R'[10,3]
			else di
			di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
		}
		
		if (("`st'"=="ex" | "`st'"=="co") & `nnt'>=0) {
			* PRINT NUMBER NEEDED TO TREAT for experimental/cohort studies (if asked) - R[6,.]
			if (`R'[6,4]==1 | `R'[6,4]==2) { 
				* CI positive or negative
				if ((`R'[6,4]==1 & `nnt'==1) | (`R'[6,4]==2 & `nnt'==0)) local lbl = "NNT Benefit (NNTB)"
				if ((`R'[6,4]==1 & `nnt'==0) | (`R'[6,4]==2 & `nnt'==1)) local lbl = "NNT Harm (NNTH)"
				di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[6,1] "`z' {c |}" /*
				*/ 	as res _col(46) %9.0g `R'[6,2] as txt " to " as res %9.0g `R'[6,3] as txt " Newcombe"
			}
			if (`R'[6,4]==3) {
				* CI negative AND positive
				di as txt "{ralign 21:NNT Benefit (NNTB)} {c |} " _c
				if ((`R'[1,1]>=0 & `nnt'==1) | (`R'[1,1]<=0 & `nnt'==0)) {
					if (`R'[1,1]!=0) di as res %9.0g `R'[6,1] "`z' {c |}" _c
					else di as txt " infinity  {c |}" _c		// RD = 0, NNT -> infinity
					di as res _col(46) %9.0g `R'[6,2] _c
				}
				if ((`R'[1,1]>0 & `nnt'==0) | (`R'[1,1]<0 & `nnt'==1)) {
					di as res _col(36) "{c |}" _col(46) %9.0g `R'[6,3] _c
				}
				di as txt " to  infinity Newcombe"
				
				di as txt "{ralign 21:NNT Harm (NNTH)} {c |} " _c
				if ((`R'[1,1]>=0 & `nnt'==1) | (`R'[1,1]<=0 & `nnt'==0)) {
					if (`R'[1,1]!=0) di as res _col(36) "{c |}" _c
					else di as txt " infinity  {c |}" _c		// RD = 0, NNT -> infinity
					di as res _col(46) %9.0g `R'[6,3] _c
				}
				if ((`R'[1,1]>0 & `nnt'==0) | (`R'[1,1]<0 & `nnt'==1)) {
					di as res %9.0g `R'[6,1] "  {c |}" _col(46) %9.0g `R'[6,2] _c
				}
				di as txt " to  infinity Newcombe"
			}
			di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
		}
		
		* PR/RR: Prevalence / Risk Ratio - R[3,.]
		di as txt _col(23) "{c |} " _col(36) "{c |}  SE(ln" cond("`st'"=="cs","PR","RR") ")"
		local lbl = cond("`st'"=="cs","Prev. Ratio (PR)","Risk Ratio (RR)")
		di as txt "{ralign 21: `lbl'} {c |} " as res %9.0g `R'[3,1] as txt "`z' {c |} " /*
		*/	as res %9.0g `R'[3,4] "  " as res %9.0g `R'[3,2] " " %9.0g `R'[3,3] " Wald"
		* PR/RR reciprocal
		if (`R'[3,1]<1) di as txt "{ralign 21:reciprocal} {c |} " as res %9.0g 1/`R'[3,1] "  {c |} " /*
		*/ 		_col(49) %9.0g 1/`R'[3,3] " " %9.0g 1/`R'[3,2]
		di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
		
		* OR: Odds Ratio - R[4,.] (CI Cornfield) + R[5,.] (CI Woolf)
		local lbl = cond("`st'"=="cs","POR","OR")
		di as txt _col(23) "{c |}" _col(36) "{c |}" _col(39) "SE(ln`lbl')" _c
		if (`ldec'==0) {
			di _col(49) as res %9.0g `R'[4,2] " " %9.0g `R'[4,3] as res " Cornfield{c `note'}"
			di as txt "{ralign 21:Odds Ratio (`lbl')} {c |} " as res %9.0g `R'[4,1] "`z' {c |} " /*
			*/		%9.0g `R'[5,4] "  " %9.0g `R'[5,2] " " %9.0g `R'[5,3] as txt " Woolf"
		}
		else {
			di _n as txt "{ralign 21:Odds Ratio (`lbl')} {c |} " as res %9.0g `R'[4,1] "`z' {c |} " /*
			*/		as res %9.0g `R'[5,4] "  " %9.0g `R'[5,2] " " %9.0g `R'[5,3] as txt " Woolf"
		}
		* OR reciprocal
		if (`R'[4,1]<1) {
			di as txt "{ralign 21:reciprocal} {c |} " as res %9.0g 1/`R'[4,1] "  {c |} " _col(49) _c
			if (`ldec'==0) di as res %9.0g 1/`R'[4,3] " " %9.0g 1/`R'[4,2] " Cornfield{c `note'}"
			else di as res %9.0g 1/`R'[5,3] " " %9.0g 1/`R'[5,2] as txt " Woolf"
			if (`ldec'==0) di as res _col(23) "{c |}" _col(36) "{c |} " _col(49) /*
			*/		%9.0g 1/`R'[5,3] " " %9.0g 1/`R'[5,2] as txt " Woolf"
		}
	}
	else {
		** CASE CONTROL studies
		* OR: Odds Ratio - R[4,.] (CI Exact) + R[9,.] (CI Cornfield) + R[5,.] (CI Woolf)
		di as txt "{ralign 21:Odds Ratio (OR)} {c |} " as res %9.0g `R'[4,1] "`z' {c |} " _c
		if (`ldec'==0) {
			di as res _col(49) %9.0g `R'[4,2] " " %9.0g `R'[4,3] as res " Exact{c `note'}"
			di as res _col(23) "{c |}" _col(36) "{c |}" _col(49) %9.0g `R'[9,2] " " %9.0g `R'[9,3] as txt " Cornfield"
			di as res _col(23) "{c |}" _col(36) "{c |} " %9.0g `R'[5,4] "  " %9.0g `R'[5,2] " " %9.0g `R'[5,3] as txt " Woolf"
		}
		else {
			di as res %9.0g `R'[5,4] "  " %9.0g `R'[5,2] " " %9.0g `R'[5,3] as txt " Woolf"
		}
		* OR reciprocal
		if (`R'[4,1]<1) {
			di as txt "{ralign 21:reciprocal} {c |} " as res %9.0g 1/`R'[4,1] "  {c |} " _col(49) _c
			if (`ldec'==0) {
				di as res %9.0g 1/`R'[4,3] " " %9.0g 1/`R'[4,2] " Exact{c `note'}"
				di as res _col(23) "{c |}" _col(36) "{c |} " _col(49) %9.0g 1/`R'[9,3] " " %9.0g 1/`R'[9,2] as txt " Cornfield"
				di as res _col(23) "{c |}" _col(36) "{c |} " _col(49) %9.0g 1/`R'[5,3] " " %9.0g 1/`R'[5,2] as txt " Woolf"
			}
			else {
				di as res %9.0g 1/`R'[5,3] " " %9.0g 1/`R'[5,2] as txt " Woolf"
			}
		}
		
		if (`r'!=0 | `pe'!=0) {
			di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
			* if r or pe provided: Pop. Risk R[1,1], PR Prop. Ratio R[2,1], PD Prop. Diff. R[3,1]
			di as txt "{ralign 21:Population Risk} {c |} " as res %9.0g `R'[1,1] "`z' {c |} " cond(`r'!=0,"(external)","(estimated)")
			di as txt "{ralign 21:Proportion Ratio (PR)} {c |} " as res %9.0g `R'[2,1] "`z' {c |}"
			di as txt "{ralign 21:Proportion Diff. (PD)} {c |} " as res %9.0g `R'[3,1] "`z' {c |}"
			
			if (`det'==1) {
				* PRINT DETAIL RESULTS (if asked)
				* PDp: Prop. Diff. pop. - R[6,1]
				di as txt "Prop. Diff. pop. (PDp){c |} " as res %9.0g `R'[6,1] "`z' {c |}"
				di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
				
				* AF/PFe: Attributable (R[7,4]=1) / Preventable (R[7,4]=2) fraction in exposed - R[7,.]
				local lbl = cond(`R'[7,4]==1,"Attr. Frac. exp.","Prev. Frac. exp.")
				di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[7,1] "`z' {c |} " as txt "Std. Err." _col(49) _c
				if (`R'[7,2]<. | `R'[7,3]<.) di as res  %9.0g `R'[7,2] " " %9.0g `R'[7,3] " " as txt cond(`ldec'==0,"Exact","Woolf")
				else di ""
			
				* AF/PFp: Attributable (R[7,4]=1) / Preventable (R[7,4]=2) fraction in population - R[8,.]
				local lbl = cond(`R'[7,4]==1,"Attr. Frac. pop.","Prev. Frac. pop.")
				di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[8,1] "`z' {c |} " _c
				if (`R'[8,4]<.) di as res %9.0g `R'[8,4] "  " _c
				else di _col(49) _c
				if (`R'[8,2]<. | `R'[8,3]<.) di as res  %9.0g `R'[8,2] " " %9.0g `R'[8,3]
				else di ""
			}
		}
	}
	
	*** PRINT ASSOCIATION 
	local lbl = cond(`chi2_type'==1,"Pearson Chi2","Mantel-Haenszel")
	local n = cond(`C'[3,1]==1,"**","")

	di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
	di as txt "Association" _col(23) "{c |}" _col(36) "{c |}"
	di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `C'[1,1] "`n'" _col(36) "{c |}" /*
	*/	as txt " p= " as res %6.4f `C'[1,2] as txt " (2-sided)"
	di as txt "{ralign 21:Corrected} {c |} " as res %9.0g `C'[2,1] _col(36) "{c |}" /*
	*/	as txt " p= " as res %6.4f `C'[2,2] as txt " (2-sided)"
	if (`ldec'==0 | (`ldec'==1 & `lzero'==1)) di as txt "{ralign 21:Fisher Exact Test} {c |} " _col(36) "{c |}" /*
	*/	as txt " p= " as res %6.4f `C'[3,2] as txt " (2-sided)"
	di as txt "{hline 22}{c BT}{hline 12}{c BT}{hline 20}"
	
	* PRINT NOTES/WARNINGS
	di as txt "{c `note'}Recommended CI"
	if (`lzero') {
		if ("`zero'"=="c") di as txt "(*)Computed with a constant continuity correction (k=0.5)"
		if ("`zero'"=="p") di as txt "(*)Computed with a proportional continuity correction"
		if ("`zero'"=="r") di as txt "(*)Computed with a reciprocal continuity correction"
	}
	if (`C'[3,1]==1) di as txt "(**)WARNING: Small samples for Association Chi-Square"
	
	
	*** STORE RESULTS
	* 2x2 table data
	return scalar a0 = `D'[1,1]
	return scalar a1 = `D'[1,2]
	return scalar m1 = `D'[1,3]
	return scalar b0 = `D'[2,1]
	return scalar b1 = `D'[2,2]
	return scalar m0 = `D'[2,3]
	return scalar n0 = `D'[3,1]
	return scalar n1 = `D'[3,2]
	return scalar n = `D'[3,3]

	* estimation results
	if ("`st'"=="cs") {
		* PD: Prevalence difference 
		return scalar pd = `R'[1,1]
		return scalar lb_pd_new = `R'[1,2]
		return scalar ub_pd_new = `R'[1,3]
		return scalar se_pd_wald = `R'[2,4]
		return scalar lb_pd_wald = `R'[2,2]
		return scalar ub_pd_wald = `R'[2,3]
		* PR: Prevalence ratio 
		return scalar pr = `R'[3,1]
		return scalar se_pr = `R'[3,4]
		return scalar lb_pr = `R'[3,2]
		return scalar ub_pr = `R'[3,3]
	}
	if ("`st'"=="co" | "`st'"=="ex") {
		* RD: Risk difference 
		return scalar rd = `R'[1,1]
		return scalar lb_rd_new = `R'[1,2]
		return scalar ub_rd_new = `R'[1,3]
		return scalar se_rd_wald = `R'[2,4]
		return scalar lb_rd_wald = `R'[2,2]
		return scalar ub_rd_wald = `R'[2,3]
		* RR: Risk ratio 
		return scalar rr = `R'[3,1]
		return scalar se_rr = `R'[3,4]
		return scalar lb_rr = `R'[3,2]
		return scalar ub_rr = `R'[3,3]
	}
	* OR: Odds ratio 
	return scalar or = `R'[4,1]
	if (`ldec'==0) {
		if ("`st'"!="cc") {
			return scalar lb_or_corn = `R'[4,2]
			return scalar ub_or_corn = `R'[4,3]
		}
		if ("`st'"=="cc") {
			return scalar lb_or_exact = `R'[4,2]
			return scalar ub_or_exact = `R'[4,3]
			return scalar lb_or_corn = `R'[9,2]
			return scalar ub_or_corn = `R'[9,3]
		}
	}
	return scalar lb_or_woolf = `R'[5,2]
	return scalar ub_or_woolf = `R'[5,3]
	return scalar se_or_woolf = `R'[5,4]
	if ("`st'"=="cc" & (`r'!=0 | `pe'!=0)) {
		* PR, PD for case-control if r or pe external
		return scalar pe = `rpext'
		return scalar r = `R'[1,1]
		return scalar pr = `R'[2,1]
		return scalar pd = `R'[3,1]
	}

	if (`det'==1) {
		* detail results
		if ("`st'"!="cc") {
			* Prop. exposed pop.
			return scalar pe_pop = `R'[7,1]
			* Prevalence/Risk diff. pop.
			if ("`st'"=="cs") {
				return scalar pdp = `R'[8,1]
				return scalar se_pdp = `R'[8,4]
				return scalar lb_pdp = `R'[8,2]
				return scalar ub_pdp = `R'[8,3]
			}
			else {
				return scalar rdp = `R'[8,1]
				return scalar se_rdp = `R'[8,4]
				return scalar lb_rdp = `R'[8,2]
				return scalar ub_rdp = `R'[8,3]
			}
			*Attr./Prev. fraction exp./pop.
			if (`R'[9,4]==1) {
				return scalar afe = `R'[9,1]
				return scalar lb_afe = `R'[9,2]
				return scalar ub_afe = `R'[9,3]
				return scalar afp = `R'[10,1]
				return scalar se_afp = `R'[10,4]
				return scalar lb_afp = `R'[10,2]
				return scalar ub_afp = `R'[10,3]
			}
			if (`R'[9,4]==2) {
				return scalar pfe = `R'[9,1]
				return scalar lb_pfe = `R'[9,2]
				return scalar ub_pfe = `R'[9,3]
				return scalar pfp = `R'[10,1]
				return scalar lb_pfp = `R'[10,2]
				return scalar ub_pfp = `R'[10,3]
			}
		}
		if ("`st'"=="cc") {
			return scalar pdp = `R'[6,1]
			if (`R'[7,4]==1) {
				return scalar afe = `R'[7,1]
				return scalar lb_afe = `R'[7,2]
				return scalar ub_afe = `R'[7,3]
				return scalar afp = `R'[8,1]
				return scalar se_afp = `R'[8,4]
				return scalar lb_afp = `R'[8,2]
				return scalar ub_afp = `R'[8,3]
			}
			if (`R'[7,4]==2) {
				return scalar pfe = `R'[7,1]
				return scalar lb_pfe = `R'[7,2]
				return scalar ub_pfe = `R'[7,3]
				return scalar pfp = `R'[8,1]
				return scalar lb_pfp = `R'[8,2]
				return scalar ub_pfp = `R'[8,3]
			}
		}
	}

	* Chi2
	return scalar chi2 = `C'[1,1]
	return scalar p = `C'[1,2]
	return scalar chi2_corr = `C'[2,1]
	return scalar p_corr = `C'[2,2]
	if (`ldec'==0 | (`ldec'==1 & `lzero'==1)) return scalar p_exact = `C'[3,2]

end


program define print_error
	args message
	display in red "`message'" 
	exit 198
end
