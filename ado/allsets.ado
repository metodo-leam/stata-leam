*! version 1.3.7  12mar2026 JM. Domenech, JB.Navarro, R. Sesma

program allsets
	version 12
	syntax varlist(fv) [if] [in] [pweight], [LINear LOGistic COX FINegray /*
	*/	fixed(varlist numeric) using(string) noreplace /*
	*/	maxvar(numlist integer min=1 max=1 >0) 		   /*
	*/  minvar(numlist integer min=1 max=1 >0)		   /*
	*/	compete(string) nohierarchical nst(string)]

	tempname nrows m r

	* check type
	local type = trim("`linear' `logistic' `cox' `finegray'")
	if ("`type'"=="") print_error "missing regression type -- specify one of linear, logistic, cox, finegray"
	if (`: word count("`type'")'>1) print_error "only one of linear, logistic, cox, finegray is needed"
	if ("`linear'"=="" & "`logistic'"=="" & "`cox'"=="" & "`finegray'"=="") {
		print_error "missing regression type -- specify one of linear, logistic, cox, finegray"
	}
	
	* check compete --finegray
	if ("`type'"=="finegray" & "`compete'"=="") print_error "compete option missing for finegray"
	if ("`type'"!="finegray" & "`compete'"!="") print_error "compete option not needed"

	* check weight
	if ("`weight'"!="" & ("`type'"=="cox" | "`type'"=="finegray")) {
		print_error "weight not necessary for cox type; include weight on the stset command"
	}

	* check filename
	if ("`using'"=="") local using "allsets_results.dta"
	capture confirm new file `"`using'"'
	if ("`replace'"!="" & _rc==602) {
		display in red "using() invalid -- the file exists and the noreplace option is specified"
		exit 198
	}
	else {
		if (_rc==603) {
			display in red "using() invalid -- the filename is invalid, the specified directory does not exist,"
			display in red "or the directory permissions do not allow to create a new file"
			exit 198
		}
	}
	if ("`hierarchical'"!="") local hierarch 0
	else local hierarch 1

	marksample touse, novarlist			// ifin marksample

	if ("`maxvar'"=="") local maxvar 0		// maximum number of variables; ALL by default
	if ("`minvar'"=="") local minvar 1		// minimum number of variables; 1 by default
	if (`maxvar'>0 & `minvar'>`maxvar') print_error "minvar option can't be greater than maxvar"

	* dependent & independent variables
	if ("`type'"!="cox" & "`type'"!="finegray") gettoken dep vars: varlist
	else local vars `varlist'
	local vars: list uniq vars		// remove possible duplicate varnames
	local error 0
	if ("`type'"!="cox" & "`type'"!="finegray") local error: list dep in vars	// dependent var can't be one of the independent
	if (`error') print_error "dependent variable `dep' can't be one of the independent variables"

	* loop to get the independent and interaction variables
	local indep ""
	local inter ""
	local cat ""
	local cont ""
	foreach var in `vars' {
		if (strpos("`var'","#")==0){		// non interaction terms
			cleanvarname `var'
			local indep = "`indep'" + " " + "`r(v_name)'"
			if (strpos("`var'",".")==0) local cont = "`cont'" + " " + "`r(v_name)'"
			else local cat = "`cat'" + " " + "`r(v_name)'"
		}
		else {		// interaction terms
			local inter = "`inter'" + " " + "`var'"
		}
	}
	* loop to verify the defined interactions
	if (`hierarch') {
		foreach in in `inter' {
			local temp : subinstr local in "#" " "
			if (strpos("`temp'","#")>0) print_error "`in' is not a valid interaction"

			local terms ""
			local nlen : word count `temp'
			forvalues i = 1/`nlen'{
				local t : word `i' of `temp'
				cleanvarname `t'
				local terms = "`terms'" + " `r(v_name)'"
			}

			local ok : list terms in indep
			if (`ok'==0) print_error "`in' is not a valid interaction, terms missing in the model"
		}
	}
	* loop to verify fixed variables
	if ("`fixed'"!="") {
		local fixed_vars
		foreach var in `fixed' {
			local is_cont: list var in cont
			local is_cat: list var in cat
			if (!`is_cont' & !`is_cat') print_error "fixed variable `var' missing in the model"
			foreach v in `vars' {
				if (strpos("`v'","#")==0){		// non interaction terms
					cleanvarname `v'
					if ("`var'"=="`r(v_name)'") local fixed_vars `fixed_vars' `v'
				}
			}
		}
	}

	* build conditions 
	local cond
	* fixed variables
	if ("`fixed'"!="") {
		foreach v in `fixed' {
			if (strpos("`v'","#")==0) {		// fixed vars must be non interaction terms
				getvarnum `v', vars(`vars')
				local cond "`cond' `r(num)';.;."
			}
		}
	}
	* interactions (hierarchical models)
	if (`hierarch') {
		foreach i in `inter' {
			local terms = subinstr("`i'","#"," ",1)
			local terms = subinstr(subinstr("`terms'","c.","",.),"i.","",.)
			local cond "`cond' "
			foreach t in `terms' {
				getvarnum `t', vars(`vars')
				local cond "`cond'`r(num)';"
			}
			getvarnum `i', vars(`vars')
			local cond "`cond'`r(num)'"
		}
	}

	if ("`type'"!="cox" & "`type'"!="finegray") local lvars = "`dep' `indep'"		// variable list
	else local lvars = "`indep'"

	quietly count if `touse'
	local nobs = r(N)					// total number of observations
	if (`nobs'==0) print_error "no observations"

	* variable list
	if ("`type'"!="cox" & "`type'"!="finegray") local lvars = "`dep' `indep'"		
	else local lvars = "`indep'"
	* count number of valid values of each var
	tempname N
	qui tabstat `lvars' if `touse', stat(count) save
	matrix `N' = r(StatTotal)

	markout `touse' `lvars'				// exclude missing values of list vars
	quietly count if `touse'
	local nvalid = r(N)					// number of valid cases

	preserve
	quietly keep if `touse'				// select cases a priori

	if ("`type'"=="cox" | "`type'"=="finegray") {
		capture which somersd			// verify somersd user-command is installed
		local somersd = (_rc==0)
		if (!`somersd') {
			di ""
			di in red "{bf:WARNING!}"
			di in red "  User-defined command {bf:somersd} is not installed."
			di in red "  This command is necessary to compute Harrell's C with Time Variable Covariates."
			if (`c(stata_version)'>=16) di in red "  Execute {bf:ssc install somersd} to install."
			else di in red `"  Execute {bf:net install somersd, force from("http://www.rogernewsonresources.org.uk/stata12")} to install."'

			exit 198
		}
	}

	* get results in Mata
	mata: getresults("`dep'", "`vars'", "`inter'","`type'",`hierarch',"`cond'","`fixed_vars'", /*
	*/	`minvar',`maxvar',"`weight'","`exp'","`compete'")
	
	* Mata code creates a new dataset with the results
	qui compress

	* properties and format
	label variable NVar "Number of variables"
	label variable Variables "Variables"
	qui replace Variables = trim(Variables)
	label variable AIC "Akaike Information Criterion"
	label variable BIC "Schwarz Bayesian Criterion"
	label variable pValue "Significance of model test"
	format pValue %5.3f
	if (`minvar'==1 & `maxvar'==0) qui drop pValue			// pValue only for minvar,maxvar calls
	format b %8.0g
	format AIC BIC %9.1f
	if ("`type'"=="linear") {
		if (`maxvar'==1) {
			label variable b "Coefficient"
			sort Cp iCat
			qui replace Cp = . if R2Adj==.
		}
		else {
			sort Cp
			qui drop b
		}
		label variable Cp "Mallows' Prediction Criterion"
		label variable R2 "R Square"
		label variable R2Adj "Adjusted R Square"
		format Cp %9.2f
		format R2Adj R2 %5.3f		
		local sumvars Cp R2Adj AIC BIC R2
	}
	if ("`type'"=="logistic") {
		if (`maxvar'==1) {
			rename b OR
			qui replace OR = exp(OR)
			label variable OR "Odds Ratio"
			sort AIC iCat
			qui replace AIC = . if BIC==.
		}
		else {
			sort AIC
			qui drop b
		}
		label variable AUC "Area Under the Curve"
		label variable Se "Sensibility (%)"
		label variable Sp "Specificity (%)"
		label variable pfitHL "Hosmer-Lemeshow goodness-of-fit"
		label variable pGof "Pearson goodness-of-fit"
		format AIC BIC %9.1f
		format AUC pfitHL pGof %5.3f
		format Se Sp %5.1f
		if ("`weight'"!="") qui drop pfitHL pGof		// p not computed on weighted regressions
	}
	if ("`type'"=="cox" | "`type'"=="finegray") {
		if (`maxvar'==1) {
			rename b HR
			qui replace HR = exp(HR)
			label variable HR "Hazard Ratio"
			sort AIC iCat
			qui replace AIC = . if BIC==.
		}
		else {
			sort AIC
			qui drop b
		}
		label variable HarrellC "Harrell's C"
		label variable _2ll "-2 Log likelihood"
		format AIC BIC %9.1f
		format _2ll %6.1f
		format HarrellC %5.3f
	}
	qui drop iCat

	local fmt : format Variables
	local fmt : subinstr local fmt "%" "%-"
	format Variables `fmt'
	quietly save `"`using'"', replace		// save results
	local path = c(filename)				// saved results path directory

	* count number of submodels estimated
	qui count if !missing(NVar)
	local nmodels = r(N)
	if (`maxvar'!=0) local nmodels = `nmodels'-1

	restore     // restore original data

	// print results
	di ""
	di ""
	local c = upper(substr("`type'",1,1)) + lower(substr("`type'",2,.))
	if ("`type'"=="cox") local c = "Cox proportional hazards"
	if ("`type'"=="finegray") local c = "Competing-risk"
	local c2 = cond("`type'"=="finegray","(Fine-Gray)","")
	di "ALLSETS - `c' regression `c2'"
	if ("`nst'"!="") di as txt "{bf:STUDY:} `nst'"

	di
	display as txt "ALL VARIABLES"
	if ("`type'"=="linear" | "`type'"=="logistic") display as text "Dependent: " as result "`dep'"
	display as text "Continuous: " as result "`cont'"
	if ("`cat'"!="") display as text "Categoricals: " as result "`cat'"
	if ("`inter'"!="") display as text "Interactions: " as result "`inter'"
	if ("`fixed'"!="") display as text "Fixed: " as result "`fixed'"

	* print the number of valid and missing values for each variable
	local names : colnames `N'
	local len = colsof(`N')
	di ""
	di ""
	di as txt "{ralign 12:Variable} {c |} {ralign 8:Valid} {ralign 8:Missing} "
	di as txt "{hline 13}{c +}{hline 18}"
	foreach i of numlist 1/`len' {
		local name = abbrev(word("`names'",`i'),12)
		local val = `N'[1,`i']
		local mis = `nobs' - `val'
		di as txt "{ralign 12:`name'} {c |} {ralign 8:`val'} {ralign 8:`mis'} "
	}
	di as txt "{hline 13}{c BT}{hline 18}"
	di as txt "Valid number of cases (listwise): " as res `nvalid'
	if (`hierarch') di as txt "Total number of " cond(`hierarch',"hierarchical","") " submodels estimated: " as res `nmodels'
	if (`minvar'!=1) di as txt "Minimum number of variables included in the combinations: " as res `minvar'
	if (`maxvar'>0) di as txt "Maximum number of variables included in the combinations: " as res `maxvar'
	dis as txt "Total time: " as res %4.1f `timer' as txt " seconds"

	// check collinearity
	_rmcoll `vars' if `touse'

	di
	di as txt "A new dataset has been created with the results. Save your data and execute:"
	di as txt "{cmd:use {c 34}`path'{c 34}, clear} to open the dataset with the results"

	* use tabstat to compute summary results
	preserve
	use `"`path'"', clear
	if ("`type'"=="linear") local svars Cp R2Adj AIC BIC R2
	if ("`type'"=="logistic" & "`wexp'"=="") local svars AIC BIC AUC Se Sp pfitHL pGof
	if ("`type'"=="logistic" & "`wexp'"!="") local svars AIC BIC AUC Se Sp
	if ("`type'"=="cox" | "`type'"=="finegray") local svars AIC BIC HarrellC _2ll

	if (`maxvar'!=0) tabstat `svars' if NVar <= `maxvar', statistics(min max)
	tabstat `svars', statistics(min max)
	if (`maxvar'!=0) di as txt "Results computed including maximum model results"
	if ("`type'"=="logistic" & "`wexp'"!="") di as txt "For weighted data the Hosmer-Lemeshow test is not computed"
	
	if (`maxvar'!=0) {
		* for maxvar calls, drop maximum model obs
		qui drop if NVar>`maxvar' & NVar<.
		qui save, replace
	}
	
	restore
end


program cleanvarname, rclass
	version 12
	args var
	if (strpos("`var'",".")>0){
		*If there's #. get the name without the #.
		local var : subinstr local var "." " ", all
		local var : word 2 of `var'
	}
	return local v_name = "`var'"
end

program getvarnum, rclass
	syntax anything(name=t), [vars(string)]

	local num 0
	local i 1
	foreach v in `vars' {
		if (strpos("`v'","#")==0) local v = substr("`v'",strpos("`v'",".")+1,.)
		if ("`t'"=="`v'") {
			local num = `i'
			break
		}
		local ++i
	}
	return local num = `num'
end

program define print_error
	args message
	display in red "`message'"
	exit 198
end


version 12
mata:
//# struct allsets_res
struct allsets_res
{
	string scalar		dep
	string rowvector	names
	real   matrix		cond
	real   scalar       rmse_max
	real   scalar		hierarchical
	real   scalar		tvc, somersd
	string scalar		type
	string scalar		weight
	string scalar		exp
	string scalar		compete
}
end

version 12
mata:
//# getresults
void getresults(string scalar dep, string scalar indep, string scalar inter, 
	string scalar type, real scalar hierarchical, string scalar cond, 
	string scalar fixed, real scalar minvar, real scalar maxvar, 
	string scalar weight, string scalar exp, string scalar compete)
{
	struct allsets_res scalar	r
	string matrix		C, combs, t, z
	string colvector	vnames
	real colvector		results, x
	real scalar			i, n, len

	r.dep = dep
	r.names = tokens(indep)
	r.type = type
	r.hierarchical = hierarchical
	r.cond = J(0,0,.)
	if (cond!="") {
		t = tokens(cond)
		for (i=1; i<=cols(t); i++) {
			z = tokens(t[1,i],";")
			if (rows(r.cond)==0) r.cond = strtoreal(select(z,z:!=";"))
			else r.cond = r.cond \ strtoreal(select(z,z:!=";"))
		}
	}
	// r.tvc = strtoreal(st_local("tvc"))
	r.somersd = strtoreal(st_local("somersd"))
	r.weight = ""
	if (weight!="") {
		r.weight = "[" + weight + exp + "]"
		r.exp = exp
	}
	r.compete = compete

	timer_clear()		// start timer
	timer_on(1)

	if (r.type == "linear") {
		// rmse of the maximum model: needed to compute Cp
		stata("regress " + dep + " " + indep + " " + r.weight,1)
		r.rmse_max = st_numscalar("e(rmse)")
	}
	
	// get all the possible combinations
	n = cols(r.names)
	if (maxvar>0 & maxvar<cols(r.names)) len = maxvar		// maximum number of variables; ALL by default
	else len = n
	combs = J(0,0,"")
	for (i=minvar; i<=len; i++) {
		// get all combinations of n elements taken i at a time
		C = allcomb(n,i,r)
		
		// build total
		if (rows(combs)==0) combs = C
		else combs = combs \ C
	}
	if (maxvar>0 & maxvar<cols(r.names)) combs = combs \ (indep)	// add maximum model for maxvar calls
	
	// compute regression for each combination
	executereg(results,vnames,combs,r,maxvar)

	// build the results dataset
	stata("clear")
	st_addobs(rows(results))
	if (r.type=="linear") {
		x = st_addvar(("float","str2000","float","float","float","float",
			"float","float","float","float","float","float"),
			("NVar","Variables","b","pValue","Cp","R2Adj","AIC","BIC","R2","iCat"))
		st_sstore(.,"Variables",vnames)
		st_store(.,"Cp",results[.,1])
		st_store(.,"R2Adj",results[.,2])
		st_store(.,"AIC",results[.,3])
		st_store(.,"BIC",results[.,4])
		st_store(.,"R2",results[.,5])
	}
	if (r.type=="logistic") {
		x = st_addvar(("float","str2000","float","float","float","float",
			"float","float","float","float","float","float"),
			("NVar","Variables","b","pValue","AIC","BIC","AUC","Se","Sp",
				"pfitHL","pGof","iCat"))
		st_sstore(.,"Variables",vnames)
		st_store(.,"AIC",results[.,1])
		st_store(.,"BIC",results[.,2])
		st_store(.,"AUC",results[.,3])
		st_store(.,"Se",results[.,4])
		st_store(.,"Sp",results[.,5])
		st_store(.,"pfitHL",results[.,6])
		st_store(.,"pGof",results[.,7])
	}
	if (r.type=="cox" | r.type=="finegray") {
		x = st_addvar(("float","str2000","float","float","float","float","float","float","float"), 
					("NVar","Variables","b","pValue","AIC","BIC","HarrellC","_2ll","iCat"))
		st_sstore(.,"Variables",vnames)
		st_store(.,"AIC",results[.,1])
		st_store(.,"BIC",results[.,2])
		st_store(.,"HarrellC",results[.,3])
		st_store(.,"_2ll",results[.,4])
	}
	st_store(.,"NVar",results[.,9])
	st_store(.,"b",results[.,10])
	st_store(.,"pValue",results[.,11])
	st_store(.,"iCat",results[.,12])

	// set macro with timer
	timer_off(1)
	st_local("timer",strofreal(timer_value(1)[1,1]))
}
end

version 12
mata:
//# allcomb
string matrix allcomb(real scalar n, real scalar k, struct allsets_res scalar r)
// Algorithm translated from Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3.  J. Gentleman.
// Author: Michael Lacy
// Date: 12 Dec 2018
{
	real matrix 	allcomb
	string matrix 	N
	real colvector	j
	real scalar		kount, nmk
	real scalar		i, L, done
	real scalar		len, ncols, nrows
	
	allcomb = J(comb(n,k), k, .) // will hold the matrix
	kount =  0
	nmk = n - k
	// Initialize j to lower limit separately, since lower limit for
	// each index depends on lower limit for previous index.
	j = J(1,k,.)
	i = 1
	j[1] = 1
	done = 0
	while (!done) {
		if (i != k) {
			for (L = i + 1;  L <= k; ++L) {
				j[L] = j[L - 1] + 1
			}
		}
		
		// Put the current combination into the matrix to be returned
		kount = kount + 1
		allcomb[kount,.] = j

		// Back to generating combinations
		// Increment the first possible index (of loop i) among indices of
		// loops k, k-1,...,1
		i = k
		while ( (j[i] >= nmk + i) & (i > 0) ) {
			i = i-1
			if (i <=0) break  // avoids evaluating j[i] when (i == 0)
		}   
		if (i > 0) j[i] = j[i] + 1
		done = (i <= 0)
	}

	// check for fixed vars
	len = rows(r.cond)
	for (i=1; i<=len; i++) {
		if (r.cond[i,2]==.) allcomb = select(allcomb,rowsum(allcomb:==r.cond[i,1]))
	}

	ncols = cols(allcomb)
	// check for hierarchical interactions
	if (r.hierarchical) {
		for (i=1; i<=len; i++) {
			if (r.cond[i,2]<.) {
				allcomb = select(allcomb,(rowsum(allcomb:!=r.cond[i,3]):==ncols :| 
					(rowsum(allcomb:==r.cond[i,3]):==1 :& 
					 rowsum(allcomb:==r.cond[i,1]):==1 :& 
					 rowsum(allcomb:==r.cond[i,2]):==1)))
			}
		}
	}
	
	// build string matrix with var names
	nrows = rows(allcomb)
	N = J(nrows,1,"")
	for (i=1; i<=nrows; i++) {
		for (j=1; j<=ncols; j++) {
			N[i] = N[i] + " " + r.names[allcomb[i,j]]
		}
	}
	return(N)
}
end

version 12
mata:
//# executereg
void executereg(real colvector results, string colvector vnames, 
	string matrix combs, struct allsets_res scalar r, real scalar maxvar)
{
	real scalar 	rss,df_m, eN, ll, ll0, cp
	real scalar		df, chi2, hc, df_r, F
	real scalar 	aic, pHL, pPearson
	real scalar     nvar
	real scalar 	i, len, ires, icat, pct, t, nrows
	string scalar 	comb, s
	string scalar	ypred, yrecodem, fr
	string scalar	hr, invhr, censind
	real matrix		freq
	real matrix		coef
	string matrix	coef_names
	real scalar 	err

	// declare results, var names matrix
	nrows = rows(combs)
	len = nrows
	if (maxvar==1) {
		for (i=1; i<=(nrows-1); i++) {
			s = combs[i]
			if (strpos(s,".")>0 & strpos(s,"#")==0) {
				// get # of unique values of categorical variable
				s = tokens(s,".")[1,3]
				len = len + rows(uniqrows(st_data(.,(s)))) - 1
			}
		}
	}
	results = J(len,12,.)
	vnames = J(len,1,"")

	ires = 1
	icat = 0
	printf("allsets running...\n")
	printf("0%%")
	displayflush()
	pct = 10
	for (i=1; i<=nrows; i++) {
		// update progress
		t = trunc(100*(ires/len))
		if (t>=(pct+10) | t==100) {
			printf("..." + strofreal(t) + "%%")
			pct = t
			displayflush()
		}

		comb = combs[i,1]
		vnames[ires] = comb
		nvar = cols(tokens(comb))
		results[ires,9] = nvar								// NVar
		if (r.type == "linear") {
			rc = _stata("regress " + r.dep + " " + comb + " " + r.weight,1)		// execute the command
			if (rc==0) {
				rss = st_numscalar("e(rss)")
				df_m = st_numscalar("e(df_m)")
				eN = st_numscalar("e(N)")
				ll = st_numscalar("e(ll)")
				df_r = st_numscalar("e(df_r)")
				F = st_numscalar("e(F)")
				cp = rss/(r.rmse_max^2) + (2*(df_m + 1) - eN)
				results[ires,1] = cp							// Cp
				results[ires,2] = st_numscalar("e(r2_a)")		// R2Adj
				results[ires,3] = -2*ll + 2*(df_m+1)			// AIC
				results[ires,4] = -2*ll + (df_m+1)*ln(eN)		// BIC
				results[ires,5] = st_numscalar("e(r2)")			// R2
				results[ires,11] = Ftail(df_m, df_r, F)			// pValue
				if (maxvar==1) {
					coef = st_matrix("e(b)")
					coef_names = st_matrixcolstripe("e(b)")
				}
			}
			else {
				// error, execute command and display error message
				rc = _stata("regress " + r.dep + " " + comb + " " + r.weight,0)
				errprintf("\nERROR!!\nmodel " + r.dep + comb + " could not be evaluated\n")
				exit(rc)
			}
		}
		
		if (r.type == "logistic") {
			rc = _stata("logit " + r.dep + " " + comb + " " + r.weight,1)	// execute the command
			if (rc==0) {
				// the command executed without error
				df_m = st_numscalar("e(df_m)")
				eN = st_numscalar("e(N)")
				ll = st_numscalar("e(ll)")
				aic = -2*ll + 2*(df_m+1)
				results[ires,1] = aic							// AIC
				results[ires,2] = -2*ll + (df_m+1)*ln(eN)		// BIC
				results[ires,11] = st_numscalar("e(p)")			// pValue
				if (maxvar==1) {
					coef = st_matrix("e(b)")
					coef_names = st_matrixcolstripe("e(b)")
				}
				
				pHL = .
				pPearson = .
				if (r.weight=="") {
					stata("lroc, nograph",1)				
					results[ires,3] = st_numscalar("r(area)")	// AUC
					stata("estat classification",1)				
					results[ires,4] = st_numscalar("r(P_p1)")	// Se
					results[ires,5] = st_numscalar("r(P_n0)")	// Sp
					if (_stata("estat gof, group(10)",1) == 0) {
						// Hosmer-Lemeshow goodness-of-fit, pfitHL
						chi2 = st_numscalar("r(chi2)")
						df = st_numscalar("r(df)")
						pHL = 1 - chi2(df,chi2)					// pfitHL
						if (nvar==1 & df==0 & pHL==.) {
							// if df = 0 and p = .,there's only 1 binary indep variable so p must be 1
							pHL = 1
						}
					}
					if (_stata("estat gof",1) == 0) {
						// Pearson goodness-of-fit, pGof
						chi2 = st_numscalar("r(chi2)")
						df = st_numscalar("r(df)")
						pPearson = 1 - chi2(df,chi2)			// pGof
						if (nvar==1 & df==0 & pPearson==.) {
							// if df = 0 and p = .,there's only 1 binary indep variable so p must be 1
							pPearson = 1
						}
					}
				}
				else {
					// weighted regression; AUC and Se, Sp are computed in a different way to apply pweight
					ypred = st_tempname()
					stata( "predict " + ypred,1)
					// execute rocreg to obtain AUC, predict variable needed
					stata( "rocreg " + r.dep + " " + ypred + " " + r.weight + ", probit ml",1)
					results[ires,3] = st_matrix("e(b)")[1,3]	// AUC (weighted)
					//Se and Sp: build the weighted table to manually compute the values
					yrecode = st_tempname()
					fr = st_tempname()
					stata("recode " + ypred + " (0/0.5=0)(0.5/1=1), gen(" + yrecode + ")",1)
					stata("tabulate " + r.dep + " " + yrecode + " [iweight" + r.exp + "], matcell(" + fr + ")",1)
					freq = st_matrix(fr)
					results[ires,4] = 100*freq[2,2]/(freq[2,1]+freq[2,2])		// Se (weighted)
					results[ires,5] = 100*freq[1,1]/(freq[1,1]+freq[1,2])		// Sp (weighted)
				}
				results[ires,6] = pHL			// pfitHL
				results[ires,7] = pPearson		// pGof				
			}
			else {
				// error, execute command and display error message
				rc = _stata("logit " + r.dep + " " + comb + " " + r.weight,0)
				errprintf("\nERROR!!\nmodel " + r.dep + comb + " could not be evaluated\n")
				exit(rc)
			}
		}
		
		if (r.type == "cox" | r.type == "finegray") {
			// execute the command
			if (r.type == "cox") rc = _stata("stcox " + comb,1)
			if (r.type == "finegray") rc = _stata("stcrreg " + comb + ", compete(" + r.compete + ")",1)
			if (rc == 0) {
				df_m = st_numscalar("e(df_m)")
				eN = st_numscalar("e(N)")
				ll = st_numscalar("e(ll)")
				chi2 = st_numscalar("e(chi2)")
				loglik = -2*ll								//AIC,BIC,loglik
				aic = loglik + 2*(df_m)
				results[ires,1] = aic						// AIC
				results[ires,2] = -2*ll + (df_m)*ln(eN)		// BIC
				results[ires,4] = -2*ll						// _2ll
				results[ires,11] = chi2tail(df_m,chi2)		// pValue
				if (maxvar==1) {
					coef = st_matrix("e(b)")
					coef_names = st_matrixcolstripe("e(b)")
				}
				// Harrell's C
				hc = .
				if (r.somersd==1) {
					// get Harrell's C with user-defined command somersd (if installed)
					hr = st_tempname()
					invhr = st_tempname()
					censind = st_tempname()
					stata("predict " + hr,1)
					stata("gen " + invhr + "=1/" + hr,1)
					stata("gen " + censind + "=1-_d if _st==1",1)
					if (_stata("somersd _t " + invhr + " if _st==1, cenind(" + censind + ") tdist transf(c)",1)==0) {
						hc = st_matrix("e(b)")[1,1]
					}
					else {
						if (st_numscalar("c(stata_version)") < 16) {
							errprintf("\nERROR!\nsomersd user-written command execution failed\n")
							errprintf("to install somersd compatible with your Stata version, execute:\n")
							errprintf("ado uninstall somersd\nnet install somersd, force from(http://www.rogernewsonresources.org.uk/stata12)\n")
							exit(3200)
						}
					}
				}
				results[ires,3] = hc	// Harrell's C
			}
			else {
				// error, execute command and display error message
				if (r.type == "cox") rc = _stata("stcox " + comb,0)
				if (r.type == "finegray") rc = _stata("stcrreg " + comb + ", compete(" + r.compete + ")",0)
				errprintf("\nERROR!!\nmodel " + comb + " could not be evaluated\n")
				exit(rc)
			}
		}
		
		if (maxvar==1 & ires<len) {
			// special case: maxvar = 1, only 1 variable models: get coefficient
			if (strpos(comb,".")>0 & strpos(comb,"#")==0) {
				// categorical variable: get coefficients of all  categories
				icat = icat + 1
				results[ires,12] = icat
				
				k = cols(coef)-1
				if (r.type == "cox" | r.type == "finegray") k = cols(coef)
				for (j=2; j<=k; j++) {
					ires = ires+1
					results[ires,10] = coef[j]
					results[ires,12] = icat + (j-1)/100
					if (r.type == "linear") results[ires,1] = cp
					if (r.type == "logistic") results[ires,1] = aic
					if (r.type == "cox" | r.type == "finegray") results[ires,1] = aic
					vnames[ires] = coef_names[j,2]
				}				
			}
			else {
				// non categorical variable
				results[ires,10] = coef[1]
			}
		}
		ires = ires+1
	}
}
end
