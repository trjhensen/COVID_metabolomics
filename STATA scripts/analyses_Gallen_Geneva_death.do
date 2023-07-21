*COVID-19 Biomarkers of death 


clear
clear mata
clear matrix 
set more off 
capture log close
set maxvar 32000

*cd G:\COVID_analyses

*generating the list of metabolites with less than 50% imputed

import excel "D:\OneDrive - National University of Ireland, Galway\Documents\stata_rerun\STATA901METS_PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED.xlsx", sheet("Sheet1") firstrow

local list_analyse=""
local n=0

foreach j of varlist carnitine-X26119{
	quietly sum `j' if AGE!=.
	*if r(N)/304>0.8
		local n=`n'+1
		local list_analyse= "`list_analyse'" + " " + "`j'"
		*
}
display `n'
clear


import excel "D:\OneDrive - National University of Ireland, Galway\Documents\stata_rerun\STATA901METS_PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED_10July2023.xlsx", sheet("Sheet1") firstrow

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

gen death=1 if Deathdue=="Yes"
replace death=0 if Deathdue=="No"
	
*Biomarkers death longitudinal - Geneva

set obs 910

local i=1
gen _var=""
gen b_coeff=.
gen CI_h=.
gen CI_l=.
gen CI_95_death=""
gen p_val_death=.
  
foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly logit death AGE sex `j' if center==2 & DaysPost==1, or
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=exp(B[1,`indB']) in `i'
	replace CI_l=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
	replace CI_h=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_death=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly logit death AGE sex `j' if center==2 & DaysPost==1, or
	quietly test `j'
	replace p_val_death=r(p) in `i'
	local i=`i'+1
	}

export excel _var CI_95_death p_val_death using "COVID_death_Geneva_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_death=""
replace p_val_death=.
local i=1


*Biomarkers death longitudinal - St Gallen


foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly logit death AGE sex `j' if center==3 & DaysPost==1, or
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=exp(B[1,`indB']) in `i'
	replace CI_l=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
	replace CI_h=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_death=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly logit death AGE sex `j' if center==3 & DaysPost==1, or
	quietly test `j'
	replace p_val_death=r(p) in `i'
	local i=`i'+1
	}

export excel _var CI_95_death p_val_death using "COVID_death_StGallen_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_death=""
replace p_val_death=.
local i=1

*Biomarkers death longitudinal - St Gallen+Geneva

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly logit death AGE sex i.center `j' if (center==3 | center==2)  & DaysPost==1, or
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=exp(B[1,`indB']) in `i'
	replace CI_l=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
	replace CI_h=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_death=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly logit death AGE sex i.center `j' if (center==3 | center==2)  & DaysPost==1, or
	quietly test `j'
	replace p_val_death=r(p) in `i'
	local i=`i'+1
	}

export excel _var CI_95_death p_val_death using "COVID_death_GallenGeneva_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_death=""
replace p_val_death=.
local i=1


* Differential trajectories dead vs. non/dead patients StGallen & geneva

gen p_val_trajectory=.
gen p_val_time=.
gen CI_95_time=""
gen p_val_global=""

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==3 | center==2) & DaysPostHospitalisation<9
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
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==3 | center==2) & DaysPostHospitalisation<9
	quietly test DaysPostHospitalisation
	replace p_val_time=r(p) in `i'
	quietly xtreg `j' AGE sex i.center death##c.DaysPostHospitalisation if (center==3 | center==2) & DaysPostHospitalisation<9
	quietly test 1.death#c.DaysPostHospitalisation
	replace p_val_trajectory=r(p) in `i'
	quietly test 1.death#c.DaysPostHospitalisation c.DaysPostHospitalisation
	*replace p_val_global=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_time p_val_time p_val_trajectory using "COVID_death_trajectory_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_time=""
replace p_val_global=.
local i=1		

*Geneva

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==2) & DaysPostHospitalisation<9
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
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==2) & DaysPostHospitalisation<9
	quietly test DaysPostHospitalisation
	replace p_val_time=r(p) in `i'
	quietly xtreg `j' AGE sex i.center death##c.DaysPostHospitalisation if (center==2) & DaysPostHospitalisation<9
	quietly test 1.death#c.DaysPostHospitalisation
	replace p_val_trajectory=r(p) in `i'
	quietly test 1.death#c.DaysPostHospitalisation c.DaysPostHospitalisation
	replace p_val_global=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_time p_val_time p_val_trajectory p_val_global using "COVID_death_trajectory_Geneva_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_time=""
replace p_val_global=.
local i=1	

*StGallen

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==3) & DaysPostHospitalisation<9
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
	quietly xtreg `j' AGE sex i.center death DaysPostHospitalisation if (center==3) & DaysPostHospitalisation<9
	quietly test DaysPostHospitalisation
	replace p_val_time=r(p) in `i'
	quietly xtreg `j' AGE sex i.center death##c.DaysPostHospitalisation if (center==3) & DaysPostHospitalisation<9
	quietly test 1.death#c.DaysPostHospitalisation
	replace p_val_trajectory=r(p) in `i'
	quietly test 1.death#c.DaysPostHospitalisation c.DaysPostHospitalisation
	replace p_val_global=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_time p_val_time p_val_trajectory p_val_global using "COVID_death_trajectory_Gallen_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_time=""
replace p_val_global=.
local i=1	

*Difference longitudinal death vs. no death StGallen+Geneva

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex i.center death if (center==3 |center==2)
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV']) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95_death=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly xtreg `j' AGE sex i.center death if (center==3 |center==2)
	quietly test death
	replace p_val_death=r(p) in `i'
	local i=`i'+1
	}

export excel _var CI_95_death p_val_death using "COVID_death_xtreg.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_death=""
replace p_val_death=.
local i=1

*Interaction Twostage effects

log using "full_twostage_ratios.log", replace
 
foreach j of varlist `list_analyse'{ 
	local a="`j'"
	foreach k of varlist `list_analyse'{ 
		local b="`k'"
		quietly xtreg `j' `k' AGE sex i.center death##c.`k' if (center==3 |center==2)  
		quietly test 1.death#c.`k'
		if r(p)<0.000001 & "`a'"!="`b'"{ 
			display "****************"
			display "`j'" 
			display "`k'"
			display r(p)
			}
		} 
	}

log close		
 
*Pairs metabolome-wide significant 


*cycloleupro<->glycolithocholatesulfate
*cholesten3one<->Noleoylserine
*sulfate<->myristoylglycerol140
*arachidoylcarnitineC20<->MU

*hydroxylysine
*Nmethyltaurine

twoway (scatter hydroxylysine Nmethyltaurine if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci hydroxylysine Nmethyltaurine if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter cycloleupro Nmethyltaurine if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci hydroxylysine glycolithocholatesulfate if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(5-hydroxylysine) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(N-methyltaurine) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction6281e_7_raw, replace)

*stearoyl2arachidonoylGPI1
*piperine

twoway (scatter stearoyl2arachidonoylGPI1 piperine if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 piperine if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter stearoyl2arachidonoylGPI1 piperine if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 piperine if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(1-stearoyl-2-arachidonoyl-GPI) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(piperine) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction6281e_7_raw, replace)

*stearoyl2arachidonoylGPI1
*octadecenedioylcarnitineC181

twoway (scatter stearoyl2arachidonoylGPI1 piperine if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(1-stearoyl-2-arachidonoyl-GPI) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(octadecenedioylcarnitine) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) 

*stearoyl2arachidonoylGPI1
*octadecenedioylcarnitineC181

twoway (scatter stearoyl2arachidonoylGPI1 piperine if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci stearoyl2arachidonoylGPI1 octadecenedioylcarnitineC181 if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(1-stearoyl-2-arachidonoyl-GPI) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(octadecenedioylcarnitine) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) 

*myristoyl2arachidonoylGPC
*sphingomyelind181200d161
twoway (scatter myristoyl2arachidonoylGPC sphingomyelind181200d161 if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci myristoyl2arachidonoylGPC sphingomyelind181200d161 if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter myristoyl2arachidonoylGPC sphingomyelind181200d161 if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci myristoyl2arachidonoylGPC sphingomyelind181200d161 if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(1-myristoyl-2-arachidonoyl-GPC) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(sphingomyelin) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction1_raw, replace)

twoway (scatter cholesten3one Noleoylserine if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci cholesten3one Noleoylserine if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter cholesten3one Noleoylserine if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci cholesten3one Noleoylserine if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(4-cholesten-3-one) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle(N-oleoylserine) xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction2_raw, replace)


twoway (scatter sulfate myristoylglycerol140 if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci sulfate myristoylglycerol140 if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter sulfate myristoylglycerol140 if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci sulfate myristoylglycerol140 if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle(Sulfate) ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle("1-myristoylglycerol (14:0)") xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction3_raw, replace)

twoway (scatter arachidoylcarnitineC20 MU if death==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci arachidoylcarnitineC20 MU if death==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter arachidoylcarnitineC20 MU if death==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci arachidoylcarnitineC20 MU if death==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==3 | center==2, ytitle("Arachidoylcarnitine (C20)") ytitle(, size(vlarge)) ylabel(, labsize(medlarge) nogrid) xtitle("2-stearoyl-GPE (18:0)") xtitle(, size(vlarge)) xlabel(, labsize(medlarge)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction4_raw, replace)

graph combine Fig_interaction1_raw.gph Fig_interaction2_raw.gph Fig_interaction3_raw.gph Fig_interaction4_raw.gph, col(2) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction_combined_raw, replace)

