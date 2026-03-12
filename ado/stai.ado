*! version 1.2.4  15sep2025 JM. Domenech, R. Sesma
/*
Association Measures - immediate data
*/

program define stai
	version 12
	syntax anything(id="argument numlist"), /*
	*/	[Data(string) ST(string) Level(numlist max=1 >50 <100) Wilson Exact WAld 		/*
	*/	 Pearson MH R(numlist max=1 >0 <1) PE(numlist max=1 >0 <1) DETail notables		/*
	*/	 NNT(numlist integer max=1 >=0 <=1) Zero(string) rare RELatsymm nst(string)]

	* check data
	if ("`data'"=="") local data "freq"
	if ("`data'"!="freq" & "`data'"!="pt" & "`data'"!="paired") print_error "data() invalid -- invalid value"
	
	* get options
	gettoken n opt: 0, parse(",")
	
	if ("`data'"=="freq") sta__freq `anything' `opt'
	if ("`data'"=="paired") sta__pair `anything' `opt'
	if ("`data'"=="pt") sta__time `anything' `opt'

end


program define print_error
	args message
	display in red "`message'" 
	exit 198
end
