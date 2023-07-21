*Analyses Ticino severe vs moderate cases


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


import excel "D:\OneDrive - National University of Ireland, Galway\Documents\stata_rerun\STATA901METS_PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED.xlsx", sheet("Sheet1") firstrow

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
gen severe_moderate=0 if study_group==1
replace severe_moderate=1 if study_group==2
	
*Biomarkers moderate severe Ticino
set obs 2000

local i=1
gen _var=""
gen b_coeff=.
gen CI_h=.
gen CI_l=.
gen CI_95_severe=""
gen p_val_severe=.


foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex severe_moderate if (center==4)
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
	replace CI_95_severe=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	quietly  xtreg `j' AGE sex severe_moderate if (center==4)
	quietly test severe_moderate
	replace p_val_severe=r(p) in `i'
	local i=`i'+1
	}

export excel _var CI_95_severe p_val_severe using "COVID_severe_ticino_xtreg.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace CI_95_severe=""
replace p_val_severe=.
local i=1

*Differential trajectory moderate-severe over time

gen p_val_trajectory=.
gen p_val_time=.
gen CI_95_time=""
gen p_val_global=.

foreach j of varlist `list_analyse'{
	local lab: variable label `j'
	replace _var="`lab'" in `i'
	quietly xtreg `j' AGE sex severe_moderate DaysPostHospitalisation if (center==4) & DaysPostHospitalisation<10
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
	quietly xtreg `j' AGE sex severe_moderate DaysPostHospitalisation if (center==4) & DaysPostHospitalisation<10
	quietly test DaysPostHospitalisation
	replace p_val_time=r(p) in `i'
	quietly xtreg `j' AGE sex severe_moderate##c.DaysPostHospitalisation if (center==4) & DaysPostHospitalisation<10 
	quietly test 1.severe_moderate#c.DaysPostHospitalisation
	replace p_val_trajectory=r(p) in `i'
	quietly test 1.severe_moderate#c.DaysPostHospitalisation c.DaysPostHospitalisation
	replace p_val_global=r(p) in `i'
	local i=`i'+1
	}
export excel _var CI_95_time p_val_time p_val_trajectory using "COVID_severe_moderate_trajectory_Ticino_re.xlsx", replace firstrow(variables)

replace _var=""
replace b_coeff=.
replace CI_h=.
replace CI_l=.
replace p_val_trajectory=.
replace p_val_time=.
replace CI_95_time=""
replace p_val_global=.
local i=1	


log using "full_twostage_ratios_moderate_severe.log", replace
 
foreach j of varlist `list_analyse'{ 
	local a="`j'"
	foreach k of varlist `list_analyse'{ 
		local b="`k'"
		quietly xtreg `j' `k' AGE sex severe_moderate##c.`k' if (center==4)  
		quietly test 1.severe_moderate#c.`k'
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

****************
*Nacetylglutamate
*cinnamoylglycine
****************
*N1methylinosine
*hydroxymyristate
****************
*N1methylinosine
*isobutyrylcarnitineC4
****************
****************
*N1methylinosine
*palmitoylGPC160
****************
*N1methylinosine
*cis4decenoate101n6
****************
*N1methylinosine
*docosahexaenoylcholine
****************
*N1methylinosine
*picolinoylglycine
****************
*5,6-dihydrouridine
*Nacetylkynurenine2
****************
*stearoylGPG180
*cinnamoylglycine
****************
*tiglylcarnitineC51DC
*Nacetylkynurenine2
****************
*cinnamoylglycine
*stearoylGPG180
****************
*sedoheptulose
*Nacetylmethionine
****************
*1-(1-enyl-palmitoyl)-GPC (P-16:0)*
*Npalmitoylsphingosined1811
****************
*(14 or 15)-methylpalmitate (a17:0 or i17:0)
*hydroxybutyrateBHBA
****************
*1-(1-enyl-stearoyl)-2-arachidonoyl-GPE (P-18:0/20:4)*
*cinnamoylglycine
****************
*formylindole
*cinnamoylglycine
****************
*allantoin
*N1methyl2pyridone5carboxami
****************
*bilirubinZZ
*androsteronesulfate
****************
*bilirubinZZ
*hydroxy2ethylpropionate
****************
*pristanate
*cinnamoylglycine
****************
*pyruvate
*thymolsulfate
****************

*Graphs
cd "G:\COVID_analyses\Results\interaction_mod_severe"

*Nacetylglutamate
*cinnamoylglycine
twoway (scatter Nacetylglutamate cinnamoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci Nacetylglutamate cinnamoylglycin if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter Nacetylglutamate cinnamoylglycin if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci Nacetylglutamate cinnamoylglycin if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N-acetyl-glutamate") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle(Cinnamoylglycine) xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction1_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction1_mod_sev_raw.png", as(png) replace
*N1methylinosine
*hydroxymyristate
twoway (scatter N1methylinosine hydroxymyristate if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine hydroxymyristate if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine hydroxymyristate if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine hydroxymyristate if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("3-hydroxymyristate") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction2_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction2_mod_sev_raw.png", as(png) replace
*N1methylinosine
*isobutyrylcarnitineC4
twoway (scatter N1methylinosine isobutyrylcarnitineC4 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine isobutyrylcarnitineC4 if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine isobutyrylcarnitineC4 if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine isobutyrylcarnitineC4 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Isobutyrylcarnitine (C4)") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction3_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction3_mod_sev_raw.png", as(png) replace
*N1methylinosine
*palmitoylGPC160
twoway (scatter N1methylinosine palmitoylGPC160 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine palmitoylGPC160 if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine palmitoylGPC160 if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine palmitoylGPC160 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("1-palmitoyl-GPC (16:0)") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction4_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction4_mod_sev_raw.png", as(png) replace
*N1methylinosine
*cis4decenoate101n6
twoway (scatter N1methylinosine cis4decenoate101n6 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine cis4decenoate101n6 if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine cis4decenoate101n6 if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine cis4decenoate101n6 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("cis-4-decenoate (10:1n6)*") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction5_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction5_mod_sev_raw.png", as(png) replace
*N1methylinosine
*docosahexaenoylcholine
twoway (scatter N1methylinosine docosahexaenoylcholine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine docosahexaenoylcholine if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine docosahexaenoylcholine if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine docosahexaenoylcholine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Docosahexaenoylcholine ") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction6_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction6_mod_sev_raw.png", as(png) replace
****************
*N1methylinosine
*picolinoylglycine
twoway (scatter N1methylinosine picolinoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine picolinoylglycine if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter N1methylinosine picolinoylglycine if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci N1methylinosine picolinoylglycine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("N1-methylinosine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Picolinoylglycine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction7_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction7_mod_sev_raw.png", as(png) replace
******* replace*********
*5,6-dihydrouridine
*Nacetylkynurenine2
twoway (scatter FQ Nacetylkynurenine2 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci FQ Nacetylkynurenine2  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter FQ Nacetylkynurenine2  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci FQ Nacetylkynurenine2 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("5,6-dihydrouridine") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("N-acetylkynurenine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction8_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction8_mod_sev_raw.png", as(png) replace
******************
*stearoylGPG180
*cinnamoylglycine
twoway (scatter stearoylGPG180 cinnamoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci stearoylGPG180 cinnamoylglycine if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter stearoylGPG180 cinnamoylglycine  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci stearoylGPG180 cinnamoylglycine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("1-stearoyl-GPG (18:0)") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Cinnamoylglycine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction9_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction9_mod_sev_raw.png", as(png) replace
****************
*tiglylcarnitineC51DC
*Nacetylkynurenine2
twoway (scatter tiglylcarnitineC51DC Nacetylkynurenine2 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci tiglylcarnitineC51DC Nacetylkynurenine2  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter tiglylcarnitineC51DC Nacetylkynurenine2  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci tiglylcarnitineC51DC Nacetylkynurenine2 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("tiglylcarnitine (C5:1-DC)") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("N-acetylkynurenine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction10_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction10_mod_sev_raw.png", as(png) replace
****************
*sedoheptulose
*Nacetylmethionine
twoway (scatter sedoheptulose Nacetylmethionine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci sedoheptulose Nacetylmethionine  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter sedoheptulose Nacetylmethionine  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci sedoheptulose Nacetylmethionine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Sedoheptulose") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("N-acetyl-methionine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction11_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction11_mod_sev_raw.png", as(png) replace
****************
*1-(1-enyl-palmitoyl)-GPC (P-16:0)*
*Npalmitoylsphingosined1811
twoway (scatter LR Npalmitoylsphingosined1811 if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci LR Npalmitoylsphingosined1811  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter LR Npalmitoylsphingosined1811  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci LR Npalmitoylsphingosined1811 if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("1-(1-enyl-palmitoyl)-GPC (P-16:0)*") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("N-palmitoyl-sphingosine (d18:1/16:0)") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction12_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction12_mod_sev_raw.png", as(png) replace
****************
*(14 or 15)-methylpalmitate (a17:0 or i17:0)
*hydroxybutyrateBHBA
twoway (scatter LW hydroxybutyrateBHBA if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci LW hydroxybutyrateBHBA   if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter LW hydroxybutyrateBHBA  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci LW hydroxybutyrateBHBA if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("(14 or 15)-methylpalmitate (a17:0 or i17:0)") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("3-hydroxybutyrate") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction13_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction13_mod_sev_raw.png", as(png) replace
****************
*1-(1-enyl-stearoyl)-2-arachidonoyl-GPE (P-18:0/20:4)*
*cinnamoylglycine
twoway (scatter RO cinnamoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci RO cinnamoylglycine  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter RO cinnamoylglycine if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci RO cinnamoylglycine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("1-(1-enyl-stearoyl)-2-arachidonoyl-GPE (P-18:0/20:4)*") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("cinnamoylglycine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction14_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction14_mod_sev_raw.png", as(png) replace
****************
*formylindole
*cinnamoylglycine
twoway (scatter formylindole cinnamoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci formylindole cinnamoylglycine  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter formylindole cinnamoylglycine if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci formylindole cinnamoylglycine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("formylindole") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("cinnamoylglycine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction15_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction15_mod_sev_raw.png", as(png) replace
****************
*allantoin
*N1methyl2pyridone5carboxami
twoway (scatter allantoin N1methyl2pyridone5carboxami if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci allantoin N1methyl2pyridone5carboxami  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter allantoin N1methyl2pyridone5carboxami if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci allantoin N1methyl2pyridone5carboxami if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Allantoin") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("N1-methyl-2-pyridone-5-carboxamide") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction16_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction16_mod_sev_raw.png", as(png) replace
****************
*bilirubinZZ
*androsteronesulfate
twoway (scatter bilirubinZZ androsteronesulfate if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci bilirubinZZ androsteronesulfate  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter bilirubinZZ androsteronesulfate if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci bilirubinZZ androsteronesulfate if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Bilirubin (Z,Z)") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Androsterone-sulfate") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction17_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction17_mod_sev_raw.png", as(png) replace
****************
*bilirubinZZ
*hydroxy2ethylpropionate
twoway (scatter bilirubinZZ hydroxy2ethylpropionate if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci bilirubinZZ hydroxy2ethylpropionate  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter bilirubinZZ hydroxy2ethylpropionate if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci bilirubinZZ hydroxy2ethylpropionate if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Bilirubin (Z,Z)") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("3-hydroxy-2-ethylpropionate") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction18_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction18_mod_sev_raw.png", as(png) replace
****************
*pristanate
*cinnamoylglycine
twoway (scatter pristanate cinnamoylglycine if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci pristanate cinnamoylglycine  if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter pristanate cinnamoylglycine if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci pristanate cinnamoylglycine if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Pristanate") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Cinnamoylglycine") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction19_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction19_mod_sev_raw.png", as(png) replace
****************
*pyruvate
*thymolsulfate
twoway (scatter pyruvate thymolsulfate if severe_moderate==1, mcolor(navy) msize(medium) msymbol(circle_hollow)) (lfitci pyruvate thymolsulfate   if severe_moderate==1, lcolor(dknavy) clwidth(thick) ciplot(rline) lcolor(dknavy) blwidth(medium) blpattern(longdash_shortdash)) (scatter pyruvate thymolsulfate  if severe_moderate==0, mcolor(maroon) msize(medium) msymbol(circle_hollow)) (lfitci pyruvate thymolsulfate  if severe_moderate==0, lcolor(cranberry) clwidth(thick) ciplot(rline) lcolor(cranberry) blpattern(longdash_shortdash)) if center==4, ytitle("Pyruvate") ytitle(, size(vlarge)) ylabel(, labsize(medsmall) nogrid) xtitle("Thymolsulfate") xtitle(, size(vlarge)) xlabel(, labsize(medsmall)) legend(off) xsize(8) ysize(8) graphregion(fcolor(white) lcolor(white)) saving(Fig_interaction20_mod_sev_raw, replace)
graph export "G:\COVID_analyses\Results\interaction_mod_severe\Fig_interaction20_mod_sev_raw.png", as(png) replace

graph combine Fig_interaction1_mod_sev_raw.gph Fig_interaction2_mod_sev_raw.gph Fig_interaction3_mod_sev_raw.gph Fig_interaction4_mod_sev_raw.gph Fig_interaction5_mod_sev_raw.gph Fig_interaction6_mod_sev_raw.gph Fig_interaction7_mod_sev_raw.gph Fig_interaction9_mod_sev_raw.gph Fig_interaction11_mod_sev_raw.gph Fig_interaction12_mod_sev_raw.gph Fig_interaction13_mod_sev_raw.gph Fig_interaction14_mod_sev_raw.gph Fig_interaction15_mod_sev_raw.gph Fig_interaction16_mod_sev_raw.gph Fig_interaction17_mod_sev_raw.gph Fig_interaction18_mod_sev_raw.gph Fig_interaction19_mod_sev_raw.gph Fig_interaction20_mod_sev_raw.gph, col(6) xsize(15) ysize(5) saving(Fig_interaction_combined_mod_sev_raw, replace) graphregion(fcolor(white) lcolor(white))  