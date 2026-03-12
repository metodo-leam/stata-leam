*! version 1.2.4  15sep2025 JM. Domenech, R. Sesma
/*
Association Measures: PERSON-TIME DATA
*/

program define sta__time, rclass
	version 12
	syntax [anything], [Data(string) Level(numlist max=1 >50 <100)  /*
	*/	 PE(numlist max=1 >0 <1) DETail notables nst(string)		/*
	*/   _incall _data(name) _valid(integer 0) _res(varname) 		/*
	*/	 _exp(varname) _time(varname) _touse(varname)]
	
	if ("`level'"=="") local level 95		//	default level: 95
	if ("`pe'"=="") local pe = 0
	local det = ("`detail'"!="")
	
	
	*** GET DATA
	* a1 a0 t1 t0
	if ("`_incall'"!="") {
		local a0 = `_data'[1,1]
		local a1 = `_data'[2,1]
		local t0 = `_data'[1,2]
		local t1 = `_data'[2,2]
	}
	else {
		tokenize `anything'
		confirm integer number `1'
		confirm integer number `2'
		confirm number `3'
		confirm number `4'
		local a1 = `1'
		local a0 = `2'
		local t1 = `3'
		local t0 = `4'
	}
	* m1 m0 n0 n1 n
	local m1 = `a0' + `a1'
	local t = `t0' + `t1'
	* data matrix
	tempname D
	matrix `D' = (`a0',`a1',`m1' \ `t0',`t1',`t')
	
	
	*** GET RESULTS
	sta__utils is_dec, data(`D')			// is there decimal values?
	local ldec = r(dec)
	sta__utils get_chi2_pt, data(`D')
	
	* get chi2 results
	tempname C
	sta__utils get_chi2_pt, data(`D')
	matrix define `C' = r(chi2)
	
	* get incidence rates
	tempname I
	sta__utils get_incidence_rates, d(`D') level(`level')
	matrix `I' = r(i)
	
	* get estimations
	tempname R
	sta__utils get_coi_results, d(`D') level(`level') pe(`pe') detail(`det')
	matrix define `R' = r(results)
	
	
	*** PRINT RESULTS
	* print title
	di as res "ASSOCIATION MEASURES: COHORT STUDY (RATE)"
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
			* get variable labels for exposure, response and time var
			sta__utils get_varlabel `_exp', len(21)
			local col_header = r(label)
			sta__utils get_varlabel `_res', len(17)
			local row_lb1 = r(label)
			sta__utils get_varlabel `_time', len(17)
			local row_lb2 = r(label)
		}
		else {
			* default labels for immediate calls
			local row_header
			local col_header
			local col_lb1 "Unexposed "
			local col_lb2 "Exposed "
			local row_lb1 "Cases"
			local row_lb2 "Person-time"
		}
		
		* header
		di
		if ("`_incall'"!="") di as txt _col(19) "{c |}{bf:{center 21:`col_header'}}" _col(41) "{c |}"
		di as txt "{bf:{ralign 17:`row_header'}} {c |}{ralign 10:`col_lb1'}{c |}{ralign 10:`col_lb2'}{c |}{ralign 9:TOTAL}"
		di as txt "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}"
		* body
		di as txt "{ralign 17:`row_lb1'} {c |} " as res %8.0g `a0' as txt " {c |} " as res %8.0g `a1' as txt " {c |} " as res %8.0g `m1'
		di as txt "{ralign 17:`row_lb2'} {c |} " as res %8.0g `t0' as txt " {c |} " as res %8.0g `t1' as txt " {c |} " as res %8.0g `t'
		di as txt "{hline 18}{c +}{hline 10}{c +}{hline 10}{c +}{hline 10}"
		* footer
		di as txt "{ralign 17:Incidence rate} {c |} " as res %8.0g `I'[1,1] as txt " {c |} " /*
		*/	as res %8.0g `I'[2,1] as txt " {c |} " as res %8.0g `I'[3,1]
 		di as txt %4.0g `level' "% CI" _col(13) "Lower {c |} " as res %8.0g `I'[1,2] as txt " {c |} " /*
		*/ 	as res %8.0g `I'[2,2] as txt " {c |} " as res %8.0g `I'[3,2]
		di as txt " (Exact)" _col(13) "Upper {c |} " as res %8.0g `I'[1,3] as txt " {c |} " /*
		*/	as res %8.0g `I'[2,3] as txt " {c |} " as res %8.0g `I'[3,3]
	}
	
	*** PRINT ESTIMATIONS
	sta__utils get_note			// get ASCII value for note Recommended CI
	local note = r(note)

	* header
	di
	di as txt "{hline 22}{c TT}{hline 12}{c TT}{hline 31}"
	di as txt "`nnt_h'{col 23}{c |}{col 26}Estimate{col 36}{c |} Std. Err.  [`level'% Conf. Interval]"
	di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"

	* ID: Incidence rate Difference - R[1,.]
	di as txt "{ralign 21:Inc. rate Diff. (ID)} {c |} " as res %9.0g `R'[1,1] as txt "  {c |} " /*
	*/	as res %9.0g `R'[1,4] "  " %9.0g `R'[1,2] " " %9.0g `R'[1,3]
	di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"

	* IR: Incidence rate Ratio - R[2,.] (CI Exact) + R[3,.] (CI Wald)
	di as txt "{ralign 21:Inc. rate Ratio (IR)} {c |} " as res %9.0g `R'[2,1] as txt "  {c |} " /*
	*/	" SE(lnIR)  " as res %9.0g `R'[2,2] " " %9.0g `R'[2,3] " Exact{c `note'}"
	di as txt _col(23) "{c |} " _col(36) "{c |} " as res %9.0g `R'[3,4] "  " %9.0g `R'[3,2] " " %9.0g `R'[3,3]

	if (`det'==1) {
		* PRINT DETAIL RESULTS (if asked)
		* AF/PFp: Attributable (R[4,4]=1) / Preventable (R[4,4]=2) fraction in exposed - R[4,.]
		local lbl = cond(`R'[4,4]==1,"Attr. Frac. exp.","Prev. Frac. exp.")
		di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[4,1] "  {c |} " _col(49) _c
		if (`R'[4,2]<. | `R'[4,3]<.) di as res  %9.0g `R'[4,2] " " %9.0g `R'[4,3] " Exact"
		di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
		* Pe: Proportion of exposed in population - R[5,1] estimated R[5,2]=1 / external R[5,2]=0
		di as txt "{ralign 21:Prop. exposed pop.} {c |} " as res %9.0g `R'[5,1] "  {c |} " /*
		*/	cond(`R'[5,2]==1,"(estimated)","(external)")
		* IDp: Inc. diff. pop. R[6,.]
		di as txt "{ralign 21:Inc. Diff. pop.(IDp)} {c |} " as res %9.0g `R'[6,1] "  {c |} " _c
		if (`R'[6,4]<.) di as res %9.0g `R'[6,4] "  " _c
		di _col(49) _c
		if (`R'[6,2]<. | `R'[6,3]<.) di as res  %9.0g `R'[6,2] " " %9.0g `R'[6,3]
		else di
		* AF/PFp: Attributable (R[4,4]=1) / Preventable (R[4,4]=2) fraction in population - R[7,.]
		local lbl = cond(`R'[4,4]==1,"Attr. Frac. pop.","Prev. Frac. pop.")
		di as txt "{ralign 21:`lbl'} {c |} " as res %9.0g `R'[7,1] "  {c |} " _c
		if (`R'[7,4]<.) di as res %9.0g `R'[7,4] "  " _c
		di _col(49) _c
		if (`R'[7,2]<. | `R'[7,3]<.) di as res  %9.0g `R'[7,2] " " %9.0g `R'[7,3]
		else di
	}


	*** PRINT ASSOCIATION 
	local n = cond(`C'[1,3]==1,"**","")
	di as txt "{hline 22}{c +}{hline 12}{c +}{hline 31}"
	di as txt "{ralign 21:Mantel-Haenszel Chi2} {c |} " as res %9.0g `C'[1,1] "`c'" _col(36) "{c |}" /*
	*/	as txt " p= " as res %6.4f `C'[1,2] as txt " (2-sided)"
	di as txt "{hline 22}{c BT}{hline 12}{c BT}{hline 10}"
	
	* PRINT NOTES/WARNINGS
	di as txt "{c `note'}Recommended CI"
	if (`C'[1,3]==1) di as txt "(**)WARNING: Small samples for Association Chi-Square"
	
	
	*** STORE RESULTS
	* 2x2 table data
	return scalar a0 = `D'[1,1]
	return scalar a1 = `D'[1,2]
	return scalar m1 = `D'[1,3]
	return scalar t0 = `D'[2,1]
	return scalar t1 = `D'[2,2]
	return scalar t = `D'[2,3]
	return scalar i0 = `I'[1,1]
	return scalar i1 = `I'[2,1]
	return scalar i = `I'[3,1]

	* ID: Incidence rate difference 
	return scalar id = `R'[1,1]
	return scalar lb_id = `R'[1,2]
	return scalar ub_id = `R'[1,3]
	* IR: Incidence rate ratio IR
	return scalar ir = `R'[2,1]
	return scalar lb_ir_exact = `R'[2,2]
	return scalar ub_ir_exact = `R'[2,3]
	return scalar se_ir = `R'[3,4]
	return scalar lb_ir = `R'[3,2]
	return scalar ub_ir = `R'[3,3]

	if (`det'==1) {
		* detail results
		* Prop. exposed pop.
		return scalar pe_pop = `R'[5,1]
		* Incidence diff. pop.
		return scalar idp = `R'[6,1]
		return scalar se_idp = `R'[6,4]
		return scalar lb_idp = `R'[6,2]
		return scalar ub_idp = `R'[6,3]
		* Attr./Prev. fraction exp./pop.
		if (`R'[4,4]==1) {
			return scalar afe = `R'[4,1]
			return scalar lb_afe = `R'[4,2]
			return scalar ub_afe = `R'[4,3]
			return scalar afp = `R'[7,1]
			return scalar se_afp = `R'[7,4]
			return scalar lb_afp = `R'[7,2]
			return scalar ub_afp = `R'[7,3]
		}
		if (`R'[4,4]==2) {
			return scalar pfe = `R'[4,1]
			return scalar lb_pfe = `R'[4,2]
			return scalar ub_pfe = `R'[4,3]
			return scalar pfp = `R'[7,1]
			return scalar lb_pfp = `R'[7,2]
			return scalar ub_pfp = `R'[7,3]
		}
	}

	* Chi2
	return scalar chi2 = `C'[1,1]
	return scalar p = `C'[1,2]
end
