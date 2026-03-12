*! version 1.2.4  15sep2025 JM. Domenech, R. Sesma
/*
Association Measures: PERSON-TIME DATA
*/

program define sta__pair, rclass
	version 12
	syntax [anything], [Data(string) Level(numlist max=1 >50 <100)  /*
	*/	 RELatsymm notables nst(string)	_incall _data(name) 		/*
	*/   _valid(integer 0) _res(varname) _exp(varname) _touse(varname)]
	
	if ("`level'"=="") local level 95		//	default level: 95	
	
	*** GET DATA
	* a0 a1 b0 b1
	if ("`_incall'"!="") {
		local a0 = `_data'[2,1]
		local b0 = `_data'[1,1]
		local a1 = `_data'[2,2]
		local b1 = `_data'[1,2]
		if ("`relatsymm'"!="") {
			local a1 = `_data'[1,2]
			local b1 = `_data'[2,2]
		}
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
		if ("`relatsymm'"!="") {
			local a1 = `3'
			local b1 = `1'
		}
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
	tempname P C R 
	sta__utils get_paired_results, d(`D') level(`level') `relatsymm'
	matrix define `P' = r(p)
	matrix define `R' = r(res)
	matrix define `C' = r(chi2)


	*** PRINT RESULTS
	* print title
	if ("`relatsymm'"!="") di as res "CONFIDENCE INTERVALS FOR MEASURES OF CHANGE (RELATIVE SYMMETRY)"
	else di as res "CONFIDENCE INTERVALS FOR MEASURES OF CHANGE (PAIRED SAMPLES)"
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
		if ("`_incall'"!="") {
			* get variable labels for exposure / response variables
			sta__utils get_var_labels `_exp', abb_name(23) abb_lbl(10)
			local col_header = r(vname)
			local col_lb1 = r(lbl0)
			local col_lb2 = r(lbl1)

			if ("`relatsymm'"=="") {
				local row_header = ""
				sta__utils get_var_labels `_res', abb_name(18) abb_lbl(10)
				local c1 = substr("`_res'",1,9)
				local c2 = substr("`_res'",10,.)
				*Swap rows 1-2 2x2 table (JMD change v1.1.7)
				local row_lb2 = r(lbl1)
				local row_lb1 = r(lbl0)
			}
			else {
				local row_header = "CHANGE of"
				local c1 = abbrev("`_res'",16)
				local c2 = abbrev("vs `_exp'",16)
				*Swap rows 1-2 2x2 table (JMD change v1.1.7)
				local row_lb2 = "Yes"
				local row_lb1 = "No"
			}
		}
		else {
			*Default labels for immediate calls
			local col_header "Response of exp. X"
			local col_lb1 "(-) "
			local col_lb2 "(+) "

			local row_header ""
			if ("`relatsymm'"=="") {
				local c1 "Response of"
				local c2 " exposure Y"
				*Swap rows 1-2 2x2 table (JMD change v1.1.7)
				local row_lb2 "(+)"
				local row_lb1 "(-)"
			}
			else {
				local c1 "CHANGE of Y"
				local c2 "   versus X"
				*Swap rows 1-2 2x2 table (JMD change v1.1.7)
				local row_lb2 "Yes"
				local row_lb1 " No"
			}
		}
		
		* header
		di
		di as txt _col(22) "{c |}{center 21: `col_header'}{c |}"
		di as txt "{lalign 20:`row_header'} {c |}{ralign 10:`col_lb1'}{c |}{ralign 10:`col_lb2'}{c |}    TOTAL  Proportion"
		di as txt "{hline 21}{c +}{hline 10}{c +}{hline 10}{c +}{hline 21}"
		* body
		if ("`relatsymm'"=="") {
			local t1 = cond("`_incall'"!="",9,11)
			local t2 = cond("`_incall'"!="",10,8)
			di as txt "{lalign `t1':`c1'} {ralign `t2':`row_lb1'} {c |} " as res %8.0g `b0' " {c |} " /*
			*/	as res %8.0g `b1' " {c |} " as res %8.0g `m0' _col(57) %9.0g `P'[1,3]
			di as txt "{lalign `t1':`c2'} {ralign `t2':`row_lb2'} {c |} " as res %8.0g `a0' " {c |} " /*
			*/	as res %8.0g `a1' " {c |} " as res %8.0g `m1' _col(57) %9.0g `P'[1,1]
		}
		else {
			di as txt "{lalign 16:`c1'} {ralign 3:`row_lb1'} {c |} " as res %8.0g `b0' " {c |} " /*
			*/	as res %8.0g `b1' " {c |} " as res %8.0g `m0'
			di as txt "{lalign 16:`c2'} {ralign 3:`row_lb2'} {c |} " as res %8.0g `a0' " {c |} " /*
			*/	as res %8.0g `a1' " {c |} " as res %8.0g `m1' _col(57) %9.0g `P'[1,1]
		}
		* footer
		di as txt "{hline 21}{c +}{hline 10}{c +}{hline 10}{c +}{hline 21}"
		di as txt "{ralign 20:TOTAL} {c |} " as res %8.0g `n0' " {c |} " as res %8.0g `n1' " {c |} " as res %8.0g `n'
		if ("`relatsymm'"=="") di as txt "{ralign 20:Proportion} {c |} " as res %8.0g `P'[1,4] " {c |} " as res %8.0g `P'[1,2] " {c |}"
		else di as txt "{ralign 20:Prop. of changes} {c |} " as res %8.0g `P'[1,2] " {c |} " as res %8.0g `P'[1,3] " {c |}"
	}
		
	*** PRINT ESTIMATIONS
	sta__utils get_note			// get ASCII value for note Recommended CI
	local note = r(note)
	
	di
	di as txt "{hline 21}{c TT}{hline 12}{c TT}{hline 32}"
	di as txt _col(22) "{c |}  Estimate  {c |} Std. Err.  [`level'% Conf. Interval]"
	di as txt "{hline 21}{c +}{hline 12}{c +}{hline 32}"
	
	if ("`relatsymm'"=="") {
		* Difference - R[1,.] (CI Exact) + R[2,.] (CI Newcombe) + R[3,.] (CI Asymp) + R[4,.] (CI Asymp. corr.)
		di as txt "{ralign 20:Difference} {c |} " as res %9.0g `R'[1,1] "  {c |}" /*
		*/	_col(48) %9.0g `R'[1,2] "  " %9.0g `R'[1,3] as txt " Exact"
		di as txt "{ralign 20:Pr(Y+)-Pr(X+)} {c |} " _col(35) "{c |}" /*
		*/	_col(48) as res %9.0g `R'[2,2] "  " %9.0g `R'[2,3] as txt " Newcombe{c `note'}"
		di as txt _col(22) "{c |}" _col(35) "{c |} " as res %9.0g `R'[3,4] /*
		*/	_col(48) %9.0g `R'[3,2] "  " %9.0g `R'[3,3] as txt " Asymptotic"
		di as txt _col(22) "{c |}" _col(35) "{c |} " as res _col(48) %9.0g `R'[4,2] "  " %9.0g `R'[4,3] as txt " A.correct."
		di as txt _col(22) "{c |}" _col(35) "{c |}"
		
		* Odds Ratio - R[5,.] (CI Exact) + R[6,.] (CI Wilson) + R[7,.] (CI Asymp)
		di as txt "{ralign 20:Odds ratio (OR)} {c |} " _c
		if (`b1'>0) di as res %9.0g `R'[5,1] "  {c |} " _col(48) %9.0g `R'[5,2] "  " %9.0g `R'[5,3] as txt " Exact"
		else di as res " infinity  {c |} " as res %9.0g `R'[5,2] "   infinity" as txt " Exact"
		di as txt "{ralign 20:(exposure-disease)} {c |}" _col(35) "{c |}  SE(lnOR)" _c
		if (`b1'>0) di as res _col(48) %9.0g `R'[6,2] "  " %9.0g `R'[6,3] as txt " Wilson{c `note'}"
		else di as res %9.0g `R'[6,2] "   infinity" di as txt " Wilson{c `note'}"
		if (`a0'>0 & `b1'>0) di as res _col(22) "{c |}" _col(35) "{c |} " %9.0g `R'[7,4] "  " /*
		*/	%9.0g `R'[7,2] "  " %9.0g `R'[7,3] as txt " Asymptotic"
		if (`R'[5,1]<1) {
			di as txt "{ralign 20:reciprocal} {c |} " as res %9.0g 1/`R'[5,1]  "  {c |} " 	/*
			*/	_col(48) %9.0g 1/`R'[5,3] "  " %9.0g 1/`R'[5,2] as txt " Exact"
			di as txt _col(22) "{c |}" _col(35) "{c |} " _col(48) as res %9.0g 1/`R'[6,3] "  " %9.0g 1/`R'[6,2] as txt " Wilson"
			di as txt _col(22) "{c |}" _col(35) "{c |} " _col(48) as res %9.0g 1/`R'[7,3] "  " %9.0g 1/`R'[7,2] as txt " Asymptotic"
		}
		if (`a0'==0) {
			di as txt "{ralign 20:reciprocal} {c |}  {bf:infinity}  {c |} " as res %9.0g 1/`R'[5,3] "   infinity" as txt " Exact"
			di as res _col(22) "{c |}" _col(35) "{c |} " %9.0g 1/`R'[6,3] "   infinity" as txt " Wilson"
		}
		di as txt "{hline 21}{c +}{hline 12}{c +}{hline 32}"
		
		* Exact simmetry
		di as txt "{ralign 20:EXACT SYMMETRY} {c |} " _col(35) "{c |}"
		di as txt "{ralign 20:McNemar matched-pair} {c |}" _col(35) "{c |} p = " as res %6.4f `C'[1,2] as txt " (Exact)"
		di as txt "{ralign 20:Chi-Square} {c |} " _c
		if (`C'[2,1]<.) di as res %9.0g `C'[2,1] as txt cond(`C'[2,3]==1,"*","") _c
		di as txt _col(35) "{c |} p = " as res %6.4f `C'[2,2]
		di as txt "{ralign 20:Corrected} {c |} " _c
		if (`C'[3,1]<.) di as res %9.0g `C'[3,1] _c
		di as txt _col(35) "{c |} p = " as res %6.4f `C'[3,2]
		di as txt "{hline 21}{c +}{hline 12}{c +}{hline 11}"
		
		* Test of Association
		di as txt "{ralign 20:TEST OF ASSOCIATION} {c |} " _col(35) "{c |}"
		di as txt "{ralign 20:OR} {c |} " _c
		if (`C'[4,1]<.) di as res %9.0g `C'[4,1] _c
		di as txt _col(35) "{c |}"
		di as txt "{ralign 20:Chi-Square} {c |} " _c
		if (`C'[5,1]<.) di as res %9.0g `C'[5,1] as txt cond(`C'[5,3]==1,"**","") _c
		di as txt _col(35) "{c |} p = " as res %6.4f `C'[5,2]
		di as txt "{ralign 20:Corrected} {c |} " _c
		if (`C'[6,1]<.) di as res %9.0g `C'[6,1] _c
		di as txt _col(35) "{c |} p = " as res %6.4f `C'[6,2]
	}
	else {
		* PD: Prop. Difference - R[1,.] (CI Wald) + R[2,.] (CI Newcombe)
		di as txt "{ralign 20:Prop. Diff. (PD)} {c |} " as res %9.0g `R'[1,1] "  {c |} " /*
		*/	as res %9.0g `R'[1,4] "  " %9.0g `R'[1,2] "  " %9.0g `R'[1,3] as txt " Wald"
		di as txt _col(22) "{c |}" _col(35) "{c |} " _col(48) /*
		*/	as res %9.0g `R'[2,2] "  " %9.0g `R'[2,3] as txt " Newcombe{c `note'}"
		
		* PR: Prop. Ratio - R[3,.]
		if (`a1'>0 & `a0'>0) di as txt _col(22) "{c |}" _col(35) "{c |}" _col(38) "SE(lnPR)"
		di as txt "{ralign 20:Prop. Ratio (PR)} {c |} " _c
		if (`a1'>0 & `a0'>0) {
			di as res %9.0g `R'[3,1] "  {c |} " %9.0g `R'[3,4] "  " /*
			*/	%9.0g `R'[3,2] "  " %9.0g `R'[3,3] as txt " Asymptotic"
			if (`R'[3,1]<1) di as txt "{ralign 20:reciprocal} {c |} " as res %9.0g 1/`R'[3,1] /*
			*/	"  {c |} " _col(48) %9.0g 1/`R'[3,3] "  " %9.0g 1/`R'[3,2]
		}	
		if (`a1'==0) di as res %9.0g 0 "  {c |}"
		if (`a1'==0) di as txt "{ralign 20:reciprocal} {c |}  infinity  {c |}"
		if (`a0'==0) di as res " infinity  {c |}"

		* OR: Odds ratio - R[4,.]
		if (`a1'>0 & `a0'>0 & `b1'>0 & `b0'>0) di as txt _col(22) "{c |}" _col(35) "{c |}" _col(38) "SE(lnOR)"
		di as txt "{ralign 20:Odds ratio (OR)} {c |} " _c
		if (`a1'>0 & `a0'>0 & `b1'>0 & `b0'>0) {
			di as res %9.0g `R'[4,1] "  {c |} " %9.0g `R'[4,4] "  " /*
			*/	%9.0g `R'[4,2] "  " %9.0g `R'[4,3] as txt " Asymptotic"
			if (`R'[4,1]<1) di as txt "{ralign 20:reciprocal} {c |} " as res %9.0g 1/`R'[4,1] /*
			*/	"  {c |} " _col(48) %9.0g 1/`R'[4,3] "  " %9.0g 1/`R'[4,2]
		}
		if (`a1'==0) di as res %9.0g 0 "  {c |}"
		if (`a1'==0) di as txt "{ralign 20:reciprocal} {c |}  infinity  {c |}"
		if (`a0'==0 | `b1'==0 | `b0'==0) di as res " infinity  {c |}"
		
		di as txt "{hline 21}{c +}{hline 12}{c +}{hline 32}"
		
		* Relative simmetry
		di as txt "{ralign 20:RELATIVE SYMMETRY} {c |} " _col(35) "{c |}"
		di as txt "{ralign 20:Chi-Square} {c |} " as res %9.0g `C'[1,1] as txt cond(`C'[1,3]==1,"**","") _c
		di as txt _col(35) "{c |} p = " as res %6.4f `C'[1,2]
		di as txt "{ralign 20:Corrected} {c |} " as res %9.0g `C'[2,1] as txt "  {c |} p = " as res %6.4f `C'[2,2]
	}
	di as txt "{hline 21}{c BT}{hline 12}{c BT}{hline 11}"
	
	* PRINT NOTES/WARNINGS
	di as txt "{c `note'}Recommended CI"
	if ("`relatsymm'"=="" & `C'[2,3]==1) di as txt "(*)WARNING: Small samples for McNemar matched-pair test"
	if (("`relatsymm'"=="" & `C'[5,3]==1) | ("`relatsymm'"!="" & `C'[1,3]==1)) /*
	*/	di as txt "(**)WARNING: Small samples for Association Chi-Square"
	
	*** STORE RESULTS
	* 2x2 table
	return scalar a10 = `D'[1,1]
	return scalar a11 = `D'[1,2]
	return scalar a00 = `D'[2,1]
	return scalar a01 = `D'[2,2]

	if ("`relatsymm'"=="") {
		* Difference
		return scalar d = `R'[1,1]
		return scalar lb_d_exact = `R'[1,2]
		return scalar ub_d_exact = `R'[1,3]
		return scalar lb_d_new = `R'[2,2]
		return scalar ub_d_new = `R'[2,3]
		return scalar lb_d_asym = `R'[3,2]
		return scalar ub_d_asym = `R'[3,3]
		return scalar se_d_asym = `R'[3,4]
		return scalar lb_d_asym_cor = `R'[4,2]
		return scalar ub_d_asym_cor = `R'[4,3]
		* OR: Odds Ratio
		return scalar or = `R'[5,1]
		return scalar lb_or_exact = `R'[5,2]
		return scalar ub_or_exact = `R'[5,3]
		return scalar lb_or_wilson = `R'[6,2]
		return scalar ub_or_wilson = `R'[6,3]
		return scalar lb_or_asym = `R'[7,2]
		return scalar ub_or_asym = `R'[7,3]
		return scalar se_or_asym = `R'[7,4]

		* Exact Simmetry
		return scalar p_McNemar = `C'[1,2]
		return scalar chi2_exact = `C'[2,1]
		return scalar p_exact = `C'[2,2]
		return scalar chi2_exact_corr = `C'[3,1]
		return scalar p_exact_corr = `C'[3,2]
		* Test of Association
		return scalar or_assoc = `C'[4,1]
		return scalar chi2_assoc = `C'[5,1]
		return scalar p_assoc = `C'[5,2]
		return scalar chi2_assoc_corr = `C'[6,1]
		return scalar p_assoc_corr = `C'[6,2]
	}
	else {
		* PD: Prop. Diff.
		return scalar pd = `R'[1,1]
		return scalar lb_pd_wald = `R'[1,2]
		return scalar ub_pd_wald = `R'[1,3]
		return scalar lb_pd_new = `R'[2,2]
		return scalar ub_pd_new = `R'[2,3]
		* PR: Prop. Ratio
		return scalar pr = `R'[3,1]
		return scalar lb_pr = `R'[3,2]
		return scalar ub_pr = `R'[3,3]
		return scalar se_pr = `R'[3,4]
		* OR: Odds Ratio
		return scalar or = `R'[4,1]
		return scalar lb_or = `R'[4,2]
		return scalar ub_or = `R'[4,3]
		return scalar se_or = `R'[4,4]
		* Relative Symmetry
		return scalar chi2 = `C'[1,1]
		return scalar p = `C'[1,2]
		return scalar chi2_cor = `C'[2,1]
		return scalar p_cor = `C'[2,2]
	}
end
