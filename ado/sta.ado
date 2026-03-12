*! version 1.2.4  15sep2025 JM. Domenech, R. Sesma
/*
Association Measures
*/

program define sta
	version 12
	syntax varlist(min=2 max=3 numeric) [if] [in], /*
	*/	[Data(string) ST(string) Level(numlist max=1 >50 <100) Wilson Exact WAld 	/*
	*/   Pearson MH R(numlist max=1 >0 <1) PE(numlist max=1 >0 <1) DETail notables	/*
	*/	 NNT(numlist integer max=1 >=0 <=1) Zero(string) rare RELatsymm nst(string)]

	* check data
	if ("`data'"=="") local data "freq"
	if ("`data'"!="freq" & "`data'"!="pt" & "`data'"!="paired") print_error "data() invalid -- invalid value"

	* get var names
	local nvars : word count `varlist'
	if (`nvars'==3 & "`data'"!="pt") print_error "wrong number of variables"
	tokenize `varlist'
	local res `1'			// response (row) variable
	local exp `2'			// exposure (column) variable
	local time `3'			// time variable (for person-time data)
	
	* check vars & data
	if ("`res'"=="`exp'") print_error "response and exposure variable must be different"
	if ("`data'"=="pt" & ("`res'"=="`time'" | "`exp'"=="`time'")) print_error "response, exposure and time variable must be different"
	if ("`res'"=="`by'" | "`exp'"=="`by'") print_error "response, exposure and stratum variable must be different"
	if ("`data'"=="pt" & "`time'"=="") print_error "time variable is needed for person-time analysis"
	if ("`data'"!="pt" & "`time'"!="") print_error "time variable is only needed for person-time analysis"

	* mark observations [if/in]
	marksample touse, novarlist
	qui count if `touse'			// count number of total cases
	local total = r(N)
	if (`total'==0) print_error "no observations"
	
	* check variable values
	qui levelsof `res', local(values)
	if ("`values'"!="0 1") print_error "response variable must binary, with values 0 1"
	qui levelsof `exp', local(values)
	if ("`values'"!="0 1") print_error "exposure variable must binary, with values 0 1"
	
	* get options
	gettoken vars opt: 0, parse(",")
	
	* get data and call subcommands to get results and print tables
	tempname d
	if ("`data'"=="freq" | "`data'"=="paired") {
		* get tabulate data
		qui tabulate `res' `exp' if `touse', matcell(`d')
		local valid = r(N)			// number of valid cases

		if ("`data'"=="freq") sta__freq `opt' _incall _data(`d') _valid(`valid') _res(`res') _exp(`exp') _touse(`touse')
		if ("`data'"=="paired") sta__pair `opt' _incall _data(`d') _valid(`valid') _res(`res') _exp(`exp') _touse(`touse')
	}
	if ("`data'"=="pt") {
		* get collapsed data
		preserve
		markout `touse' `exp' `res' `time'
		qui count if `touse'		// number of valid cases
		local valid = r(N)
		qui collapse (sum) __a=`res' __t=`time' if `touse', by(`exp')
		mkmat __a __t, matrix(`d')
		restore
		
		sta__time `opt' _incall _data(`d') _valid(`valid') _res(`res') _exp(`exp') _time(`time') _touse(`touse')
	}
	
end


program define print_error
	args message
	display in red "`message'" 
	exit 198
end
