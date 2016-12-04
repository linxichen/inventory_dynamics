
gen recession_dum = 0
replace recession_dum = 1 if  -45 <= date & date <=-41
replace recession_dum = 1 if  -26 <= date & date <=-23
replace recession_dum = 1 if  -10 <= date & date <=-7
replace recession_dum = 1 if  1 <= date & date <= 4
replace recession_dum = 1 if  39 <= date & date <= 43
replace recession_dum = 1 if  55 <= date & date <= 60
replace recession_dum = 1 if  80 <= date & date <= 82
replace recession_dum = 1 if  86 <= date & date <= 91
replace recession_dum = 1 if  122 <= date & date <= 124
replace recession_dum = 1 if  164 <= date & date <= 167
replace recession_dum = 1 if  191 <= date & date <= 197

gen regimes = 0
replace regimes = 1 if recession_dum == 1 & date < yq(1984,1)
replace regimes = 2 if recession_dum == 0 & date < yq(1984,1)
replace regimes = 3 if recession_dum == 1 & date >= yq(1984,1)
replace regimes = 4 if recession_dum == 0 & date >= yq(1984,1)

gen simple_regimes = recession_dum + 1

save "../data/quarterly_FRED.dta", replace

export delimited using "../data/quarterly_FRED.txt", datafmt replace

