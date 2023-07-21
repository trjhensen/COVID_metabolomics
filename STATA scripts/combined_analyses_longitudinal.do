*COVID-19 Data analyses longitudinal Metabolome Data (all centers together)
*Johannes Hertel - johannes.hertel@med.uni-greifswald.de

clear
clear mata
clear matrix 
set more off 
capture log close
set maxvar 32000

cd G:\COVID_analyses

*generating the list of metabolites with less than 50% imputed

import excel "G:\COVID_analyses\PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_RAW_NO_IMP.xlsx", sheet("Sheet1") firstrow

local list_analyse=""
local n=0

foreach j of varlist carnitine-X26119{
	quietly sum `j' if AGE!=.
	if r(N)/304>0.8{
		local n=`n'+1
		local list_analyse= "`list_analyse'" + " " + "`j'"
		}
	}
display `n'
clear


import excel "G:\COVID_analyses\PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED.xlsx", sheet("Sheet1") firstrow

gen s=strpos(Patient_ID, "_")-1
replace s=0 if s<0
gen id=substr(Patient_ID,1,s)
replace id=Patient_ID if s==0
order Patient_ID id
egen id_nr=group(id)
tab id_nr
gen study_group=0 if COVID_SEVERITY=="control"
replace study_group=1 if substr(COVID_SEVERITY,1,14)=="COVID_moderate"
replace study_group=2 if substr(COVID_SEVERITY,1,12)=="COVID_severe"
replace study_group=3 if substr(COVID_SEVERITY,1,10)=="Long COVID"
label define group 0 "Control" 1 "COVID Moderate" 2 "COVID severe" 3 "Long Covid"
label values study_group group
egen sex=group(GENDER)
replace sex=. if sex==3
label define gen 1 "Female" 2 "Male"
label values sex gen

gen time=strreverse(substr(strreverse(COVID_SEVERITY),1,2))
replace time="T1" if substr(time,1,1)!="T"
replace time=substr(time,2,1)
destring(time), replace
egen center=group(Location)

xtset id_nr DaysPostHospitalisation

foreach j of varlist carnitine-X26119{
	replace `j'=ln(`j')
	}

*Longitudinal association study - global differences between controls, severe and moderate cases (random effects)
cd G:\COVID_analyses\Results
gen std_gr2=study_group if study_group!=3
gen severe_moderate=0 if study_group==1
replace severe_moderate=1 if study_group==2
set obs 2000

local i=1
gen _var=""
gen b_coeff=.
gen CI_h=.
gen CI_l=.
gen CI_95_severe=""
gen CI_95_moderate=""
gen CI_95_mod_sev=""
gen p_val_global=.
gen p_val_severe=.
gen p_val_mod_sev=.
gen p_val_moderate=.

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center i.std_gr2
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_severe=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	local indV=colsof(V)-2
	local indB=colsof(B)-2
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_moderate=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center i.std_gr2
	quietly test 1.std_gr2 2.std_gr2
	replace p_val_global=r(p) in `i'
	quietly test 1.std_gr2 
	replace p_val_moderate=r(p) in `i'
	quietly test 2.std_gr2
	replace p_val_severe=r(p) in `i'
	quietly xtreg `j' AGE sex i.center severe_moderate
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_mod_sev=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center severe_moderate
	quietly test severe_moderate
	replace p_val_mod_sev=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_moderate CI_95_severe CI_95_mod_sev p_val_moderate p_val_severe p_val_mod_sev p_val_global using "COVID_study_groups_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_severe=""
replace CI_95_moderate=""
replace CI_95_mod_sev=""
replace p_val_global=.
replace p_val_severe=.
replace p_val_moderate=.
replace p_val_mod_sev=.
local i=1	
	
* Differential trajectories severe vs moderate COVID-cases 

gen p_val_trajectory=.
gen p_val_time=.
gen CI_95_time=""

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center severe_moderate DaysPostHospitalisation if DaysPostHospitalisation<9
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_time=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center severe_moderate DaysPostHospitalisation if DaysPostHospitalisation<9
	quietly test DaysPostHospitalisation
	replace p_val_time=r(p) in `i'
	quietly xtreg `j' AGE sex i.center severe_moderate##c.DaysPostHospitalisation if DaysPostHospitalisation<9
	quietly test 1.severe_moderate#c.DaysPostHospitalisation
	replace p_val_trajectory=r(p) in `i'
	quietly test 1.severe_moderate#c.DaysPostHospitalisation c.DaysPostHospitalisation
	replace p_val_global=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_time p_val_time p_val_trajectory p_val_global using "COVID_mod_sev_trajectory_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_time=""
replace p_val_global=.
local i=1		
	
*Differences severe vs moderate cases at day one of hospitalisation

gen p_val_baseline=.

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly reg `j' AGE sex i.center severe_moderate if DaysPostHospitalisation==1
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_mod_sev=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly reg `j' AGE sex i.center severe_moderate if DaysPostHospitalisation==1
	quietly test severe_moderate
	replace p_val_baseline=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_mod_sev p_val_baseline using "COVID_mod_sev_baseline.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_mod_sev=""
replace p_val_baseline=.
local i=1			
	
*Effects of location
egen center=group(Location)
gen p_val_location=.
gen CI_95_Geneva=""
gen CI_95_StGallen=""
gen CI_95_Ticino=""


foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center if study_group==2
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_Ticino=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center if study_group==2
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-2
	matrix B=e(b)
	local indB=colsof(B)-2
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_StGallen=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center if study_group==2
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-3
	matrix B=e(b)
	local indB=colsof(B)-3
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_Geneva=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center if study_group==2
	quietly test 2.center 3.center 4.center
	replace p_val_location=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_Geneva CI_95_StGallen CI_95_Ticino p_val_location using "COVID_severe_location_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_location=.
replace CI_95_Geneva=""
replace CI_95_StGallen=""
replace CI_95_Ticino=""
local i=1			
	