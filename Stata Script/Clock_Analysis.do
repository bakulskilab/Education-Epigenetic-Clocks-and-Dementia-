% Analysis of Epigenetic Clocks and Dementia 
% Cesar Higgins MD, MPH, MS; Erin B. Ware MS, PhD; Kelly M. Bakulski PhD. 


/*Import Complete Data
This sample of 4,018 participants includes cell type composition and all covariates */

use "/Users/cesarhiggins/Documents/Epigenetic Clocks /Tables/Epigenetic Paper/Code and Data/clock_complete_data.dta", clear

recode raceth (1=1) (99=.), gen(raceth_f)

label define raceth_f 1"Other/NA" 2 "Black NH" 3 "Hispanic" 4 "White NH"
label values raceth_f raceth_f 
label var raceth_f "Race"
tab raceth_f

rename PMONO monocytes
rename PLYMP lynphocytes 
gen granulocytes = PBASO+PEOS+PNEUT /*N = 33 missing values, may be people without granulocytes */

gen totalcell = lynphocytes+monocytes+granulocytes

summarize granulocytes monocytes lynphocytes totalcell

encode sex, gen(sex2) 

summarize page, detail 

encode dementia_lw, gen(dementia_lw_2)


/*Generate missing dummy variable*/

gen missing_2=.
replace missing_2=0 if dementia_lw_2!=. & edubin_low!=. & raceth_f!=. & sex2!=. & page!=. & granulocytes!=. & monocytes!=. & lynphocytes!=. & vbsi16wgtra!=. &  ///
horvath_dnamage!=. & hannum_dnamage!=. & levine_dnamage!=. & horvathskin_dnamage!=. & dnamgrimage!=. & mpoa!=. 
replace missing_2=1 if dementia_lw_2==. | edubin_low==. | raceth_f==. | sex2==. | page==. | granulocytes==. | monocytes==. | lynphocytes==. | vbsi16wgtra==. | ///
horvath_dnamage==. | hannum_dnamage==.  | levine_dnamage==. | horvathskin_dnamage==. | dnamgrimage==. | mpoa==. 
tab missing_2 

drop if missing_2==1
drop if raceth_f==1


***************************
*** CALCULATE RESIDUALS ***
***************************

* regression of each clock on chronological age 


* HORVATH 
regress horvath_dnamage page
predict horthvath_xb, xb
gen horthvath_res = horvath_dnamage-horthvath_xb
summarize horthvath_res

* residuals versus fitted values 
scatter horthvath_res horthvath_xb

* residuals versus chronological age 
scatter horthvath_res page
histogram horthvath_res , normal 

* HANNUN

regress hannum_dnamage page
predict hannum_xb, xb
gen hannum_res = hannum_dnamage-hannum_xb
summarize hannum_res

* residuals versus fitted values 
scatter hannum_res hannum_xb

* residuals versus chronological age 
scatter hannum_res page

histogram hannum_res , normal 

* LEVINE

regress levine_dnamage page
predict levine_xb, xb
gen levine_res = levine_dnamage-levine_xb
summarize levine_res

* residuals versus fitted values 
scatter levine_res levine_xb

* residuals versus chronological age 
scatter levine_res page

histogram levine_res, normal 


* HORVATH SKIN 

regress horvathskin_dnamage page
predict horvathskin_xb, xb
gen horvathskin_res = horvathskin_dnamage-horvathskin_xb
summarize horvathskin_res

* residuals versus fitted values 
scatter horvathskin_res horvathskin_xb

* residuals versus chronological age 
scatter horvathskin_res page

histogram horvathskin_res, normal 

* GRIM AGE

regress dnamgrimage page
predict grimage_xb, xb
gen grimage_res = dnamgrimage-grimage_xb
summarize grimage_res

* residuals versus fitted values 
scatter grimage_res grimage_xb

* residuals versus chronological age 
scatter grimage_res page

histogram grimage_res, normal 

* MPOA


gen mpoa_dnamage = mpoa*page 
summarize mpoa_dnamage

* MPOA

regress mpoa_dnamage page
predict mpoa_xb, xb
gen mpoa_res = mpoa_dnamage-mpoa_xb
summarize mpoa_res

* residuals versus fitted values 
scatter mpoa_res mpoa_xb

* residuals versus chronological age 
scatter mpoa_res page

histogram mpoa_res, normal 


**********************
**** NEW ANALYSIS ****
**********************


**** TABLE 1 WITH RESIDUALS AND LANGA WEIR AS THE TOP VARIABLE *** 

gen ID = _n 

svyset ID [pweight=vbsi16wgtra], strata(stratum) singleunit(certainty)  


recode dementia_lw_2 (1=1) (3=0) (2=2), gen(cog_status)
tab cog_status
la def cog_status 0 "Normal" 1 "CIND" 2 "Dementia" 
la val cog_status cog_status


/*Crosstabulation - table 1*/

tabout raceth_f sex2 edubin_low cog_status using table_1_9-27-22.xls, ///
c(col ci) f(1 1) clab(Col_% 95%_CI) svy stats(chi2) ///
npos(lab) percent ///
rep 

foreach var in page horthvath_res hannum_res levine_res horvathskin_res grimage_res mpoa_res granulocytes monocytes lynphocytes {
tabout cog_status using table_1_9-27-22.xls, append ///
c(mean `var' ci) f(1 1) clab(`var' 95%_CI) ///
sum svy stats(ttest) npos(lab) 
}

tabout powers_d cog_status using table_1.1_9-27-22.xls, ///
c(freq row col) f(0c 1p 1p) lay(cb) h1(nil) h3(nil) npos(row) ///
noffset(3) replace

***************************************
* TABLE 2 RESIDUALS DICHOTOMIZED ******
***************************************


** DICHOTOMIZING EACH CLOCK AT > TO THE MEAN 

local varlist grimage_res mpoa_res levine_res horthvath_res hannum_res  horvathskin_res 
foreach v in `varlist' {
sum `v', detail 
foreach n of numlist `r(mean)'	{
gen byte dicho_`v' = `v' > `r(mean)'
}
}
tab1 dicho_*


** Dementia models for each year of age acceleration 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit dn_lw `xvar' i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	local m `++m'
}
	esttab model_* using model_clock_res.rtf, replace ci obslast b(%9.2f) t(%9.1f) wide eform ///	
	label nonumber title("Models for Dementia") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear


** CIND models for each year of each acceleration 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit cn_lw `xvar' i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	local m `++m'
}
	esttab model_* using model_clock_cind_res.rtf, replace ci obslast b(%9.2f) t(%9.1f) wide eform ///	
	label nonumber title("Models for CIND") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear

** Dementia models for each year of age acceleration POWERS 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit powers_d `xvar' i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	local m `++m'
}
	esttab model_* using model_clock_res_powers.rtf, replace ci obslast b(%9.2f) t(%9.1f) wide eform ///	
	label nonumber title("Models for Dementia Powers") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear
	
***************************************
* TABLE 3 EDU + RESIDUAL INTERACTION **
***************************************



** AT THE MEAN LANGA-WEIRD 
** Dementia Additive Models clocks dichotomized at the mean of accelerated residuals 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit dn_lw i.edubin_low#`xvar' ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.1f)
	local m `++m'
}
	esttab model_* using model_clock_res_int_2.rtf, replace ci obslast b(%9.1f) t(%9.2f) wide eform ///	
	label nonumber title("Models for Dementia") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear


** CIND Additive Models clocks at mean of accelerated residuals 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit cn_lw i.edubin_low#`xvar' ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	local m `++m'
}
	esttab model_* using model_clock_cind_res_int_2.rtf, replace ci obslast b(%9.1f) t(%9.2f) wide eform ///	
	label nonumber title("Models for CIND") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear
	
	
** Dementia models clocks at mean of accelerated residuals  interaction 
local m=1 
	foreach xvar in dicho_grimage_res dicho_mpoa_res dicho_levine_res dicho_horthvath_res dicho_hannum_res dicho_horvathskin_res {
	eststo model_`m': svy: logit powers_d i.edubin_low#`xvar' ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	local m `++m'
}
	esttab model_* using model_clock_powers_res_int_2.rtf, replace ci obslast b(%9.1f) t(%9.2f) wide eform ///	
	label nonumber title("Models for Dementia Powers") mtitles("Grim Age" "MPOA" "Levine" "Horvath" "Hannum" "Horvath Skin")
	eststo clear
		
	


***************************************
* TABLE 4 MEDIATION ANALYSIS **********
***************************************

** INTERACTION OF GRIM CLOCK LANGA-WEIRD DEMENTIA

* create 4 levels variables 

gen edu_grim=. 
replace edu_grim=1 if dicho_grimage_res==0 & edubin_low==0
replace edu_grim=2 if dicho_grimage_res==1 & edubin_low==0
replace edu_grim=3 if dicho_grimage_res==0 & edubin_low==1
replace edu_grim=4 if dicho_grimage_res==1 & edubin_low==1
tab edu_grim

xi: svy: logit dn_lw i.edu_grim ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.1f)

	nlcom [exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1, cformat(%3.1f) /*REIR*/

	nlcom exp(_b[_Iedu_grim_4]*(-1)) - exp(_b[_Iedu_grim_3] - _b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_2] - _b[_Iedu_grim_4]) + 1, cformat(%3.1f) /*AP*/

	nlcom [exp(_b[_Iedu_grim_4]) - 1]/[exp(_b[_Iedu_grim_3]) + exp(_b[_Iedu_grim_2]) - 2], cformat(%3.1f) /*SI*/

	nlcom ln([exp(_b[_Iedu_grim_4]) - 1]) - ln([exp(_b[_Iedu_grim_3]) + exp(_b[_Iedu_grim_2]) - 2]), cformat(%3.1f) /*SI*/
	
	
## Calculating the model with survey weights for while male 


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_grimage_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=4 page=67.87 sex2=2)
	
	gen pm=.4372842
	
	gen med=.6076222-.4372842

xi: svy: logit dn_lw i.edu_grim ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_grim_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_grim_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med), cformat(%9.1f)
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_grim_3]) - 1) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom ((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))] + [((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)) + ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)
	
	

	
	
## Calculating the model with survey weights for Black male 


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_grimage_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=2 page=67.87 sex2=2)
	
	drop pm med
	
	gen pm=.6784756
	
	gen med=.8085028-.6784756

xi: svy: logit dn_lw i.edu_grim ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_grim_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_grim_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med), cformat(%9.1f)
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_grim_3]) - 1) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom ((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))] + [((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)) + ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)
	

	
## Calculating the model with survey weights for Hispanic male 


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_grimage_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=3 page=67.87 sex2=2)
	
	drop pm med
	
	gen pm=.4462461
	
	gen med=.6162456-.4462461

xi: svy: logit dn_lw i.edu_grim ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_grim_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_grim_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med), cformat(%9.1f)
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_grim_3]) - 1) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom ((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))] + [((exp(_b[_Iedu_grim_2]) - 1)*med) / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med)) + ([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med / ((exp(_b[_Iedu_grim_3]) - 1) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*pm) + (([exp(_b[_Iedu_grim_4]) - exp(_b[_Iedu_grim_3]) - exp(_b[_Iedu_grim_2])] + 1)*med) + ((exp(_b[_Iedu_grim_2]) - 1)*med))], cformat(%9.3f)	
	

******************************************************
** INTERACTION OF MPOA CLOCK LANGA-WEIRD DEMENTIA ****
******************************************************

* create 4 levels variables 

gen edu_mpoa=. 
replace edu_mpoa=1 if dicho_mpoa_res==0 & edubin_low==0
replace edu_mpoa=2 if dicho_mpoa_res==1 & edubin_low==0
replace edu_mpoa=3 if dicho_mpoa_res==0 & edubin_low==1
replace edu_mpoa=4 if dicho_mpoa_res==1 & edubin_low==1
tab edu_mpoa

xi: svy: logit dn_lw i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.1f)

	nlcom [exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1, cformat(%3.1f) /*REIR*/

	nlcom exp(_b[_Iedu_mpoa_4]*(-1)) - exp(_b[_Iedu_mpoa_3] - _b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_2] - _b[_Iedu_mpoa_4]) + 1, cformat(%3.1f) /*AP*/

	nlcom [exp(_b[_Iedu_mpoa_4]) - 1]/[exp(_b[_Iedu_mpoa_3]) + exp(_b[_Iedu_mpoa_2]) - 2], cformat(%3.1f) /*SI*/

	nlcom ln([exp(_b[_Iedu_mpoa_4]) - 1]) - ln([exp(_b[_Iedu_mpoa_3]) + exp(_b[_Iedu_mpoa_2]) - 2]), cformat(%3.1f) /*SI*/
	
	
***************************************************************
** MEDIATION DECOMPOSITION MPOA CLOCK LANGA-WEIRD DEMENTIA ****
***************************************************************
	
## Calculating the model with survey weights for while male MPOA


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=4 page=67.87 sex2=2)
	
	drop pm med
	
	gen pm=.374561
	
	gen med=.5000532-.374561

xi: svy: logit dn_lw i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	
	

	
	
## Calculating the model with survey weights for Black male 


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=2 page=67.87 sex2=2)
	
	drop pm med
	
	gen pm=.6538116 
	
	gen med=.7594729-.6538116 

xi: svy: logit dn_lw i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

	
## Calculating the model with survey weights for Hispanic male 


	svy: mean page if dn_lw!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if dn_lw!=., or
	
	margins edubin_low, at(raceth=3 page=67.87 sex2=2)
	
	drop pm med
	
	gen pm=.4595149
	
	gen med=.5861651-.4595149

xi: svy: logit dn_lw i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

***************************************************************
********** INTERACTION OF  MPOA CLOCK POWERS ******************
***************************************************************	
	
xi: svy: logit powers_d i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.1f)
	
	
	nlcom [exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1, cformat(%3.1f) /*REIR*/

	nlcom exp(_b[_Iedu_mpoa_4]*(-1)) - exp(_b[_Iedu_mpoa_3] - _b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_2] - _b[_Iedu_mpoa_4]) + 1, cformat(%3.1f) /*AP*/

	nlcom [exp(_b[_Iedu_mpoa_4]) - 1]/[exp(_b[_Iedu_mpoa_3]) + exp(_b[_Iedu_mpoa_2]) - 2], cformat(%3.1f) /*SI*/

	nlcom ln([exp(_b[_Iedu_mpoa_4]) - 1]) - ln([exp(_b[_Iedu_mpoa_3]) + exp(_b[_Iedu_mpoa_2]) - 2]), cformat(%3.1f) /*SI*/
	
	
	
***************************************************************
** MEDIATION DECOMPOSITION MPOA CLOCK POWERS ******************
***************************************************************

## Calculating the model with survey weights for while male MPOA


	svy: mean page if powers_d!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if powers_d!=., or
	
	margins edubin_low, at(raceth=4 page=67.6 sex2=2)
	
	drop pm med
	
	gen pm=.3878331
	
	gen med=.5063371-.3878331

xi: svy: logit powers_d i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	



## Calculating the model with survey weights for black male MPOA


	svy: mean page if powers_d!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if powers_d!=., or
	
	margins edubin_low, at(raceth=2 page=67.6 sex2=2)
	
	drop pm med
	
	gen pm=.6466761
	
	gen med=.7478932-.6466761

xi: svy: logit powers_d i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	



## Calculating the model with survey weights for Hispanic male MPOA


	svy: mean page if powers_d!=.
	
	svy: logit dicho_mpoa_res i.edubin_low ib(4).raceth page i.sex2 granulocytes monocytes if powers_d!=., or
	
	margins edubin_low, at(raceth=3 page=67.6 sex2=2)
	
	drop pm med
	
	gen pm=.4499018
	
	gen med=.5692968-.4499018

xi: svy: logit powers_d i.edu_mpoa ib(4).raceth page i.sex2 granulocytes monocytes, or cformat(%3.2f)
	
## CDE 

	nlcom exp(_b[_Iedu_mpoa_3]) - 1, cformat(%9.1f)	
	
## INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm, cformat(%9.1f)
	
## INT mediation
	
	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med, cformat(%9.1f)
	
## PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med, cformat(%9.1f)
	
## TOTAL 


	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med), cformat(%9.1f)
	
	
		
# Calculating the proportion mediated 


## PA-CDE 

	nlcom (exp(_b[_Iedu_mpoa_3]) - 1) / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-INT reference 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)

## PA-INT med 

	nlcom ([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## PA-PIE 

	nlcom (exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med)), cformat(%9.3f)
	
## Portion-Mediated 

	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [(exp(_b[_Iedu_mpoa_2]) - 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	

## Portion-Interaction 


	nlcom [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))] + [([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med / ((exp(_b[_Iedu_mpoa_3]) - 1) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*pm) + (([exp(_b[_Iedu_mpoa_4]) - exp(_b[_Iedu_mpoa_3]) - exp(_b[_Iedu_mpoa_2])] + 1)*med) + ((exp(_b[_Iedu_mpoa_2]) - 1)*med))], cformat(%9.3f)
	



**************************************************
************* INCLUSION/EXCLUSION ****************
**************************************************

use "/Users/cesarhiggins/Documents/Epigenetic Clocks /Tables/Epigenetic Paper/Code and Data/clock_complete_data.dta", clear

recode raceth (1=.) (99=.) , gen(raceth_f)

label define raceth_f 1"Other/NA" 2 "Black NH" 3 "Hispanic" 4 "White NH"
label values raceth_f raceth_f 
label var raceth_f "Race"
tab raceth_f

rename PMONO monocytes
rename PLYMP lynphocytes 
gen granulocytes = PBASO+PEOS+PNEUT /*N = 33 missing values, may be people without granulocytes */

gen totalcell = lynphocytes+monocytes+granulocytes

summarize granulocytes monocytes lynphocytes totalcell

encode sex, gen(sex2) 

summarize page, detail 

encode dementia_lw, gen(dementia_lw_2)



/*Generate missing dummy variable*/

gen missing_2=.
replace missing_2=0 if dementia_lw_2!=. & edubin_low!=. & raceth_f!=. & sex2!=. & page!=. & granulocytes!=. & monocytes!=. & lynphocytes!=. & vbsi16wgtra!=. &  ///
horvath_dnamage!=. & hannum_dnamage!=. & levine_dnamage!=. & horvathskin_dnamage!=. & dnamgrimage!=. & mpoa!=. 
replace missing_2=1 if dementia_lw_2==. | edubin_low==. | raceth_f==. | sex2==. | page==. | granulocytes==. | monocytes==. | lynphocytes==. | vbsi16wgtra==. | ///
horvath_dnamage==. | hannum_dnamage==.  | levine_dnamage==. | horvathskin_dnamage==. | dnamgrimage==. | mpoa==. 
tab missing_2 




* regression of each clock on chronological age to build residuals 

***************************
*** CALCULATE RESIDUALS ***
***************************

* regression of each clock on chronological age 


* HORVATH 
regress horvath_dnamage page
predict horthvath_xb, xb
gen horthvath_res = horvath_dnamage-horthvath_xb
summarize horthvath_res

* residuals versus fitted values 
scatter horthvath_res horthvath_xb

* residuals versus chronological age 
scatter horthvath_res page
histogram horthvath_res , normal 

* HANNUN

regress hannum_dnamage page
predict hannum_xb, xb
gen hannum_res = hannum_dnamage-hannum_xb
summarize hannum_res

* residuals versus fitted values 
scatter hannum_res hannum_xb

* residuals versus chronological age 
scatter hannum_res page

histogram hannum_res , normal 

* LEVINE

regress levine_dnamage page
predict levine_xb, xb
gen levine_res = levine_dnamage-levine_xb
summarize levine_res

* residuals versus fitted values 
scatter levine_res levine_xb

* residuals versus chronological age 
scatter levine_res page

histogram levine_res, normal 


* HORVATH SKIN 

regress horvathskin_dnamage page
predict horvathskin_xb, xb
gen horvathskin_res = horvathskin_dnamage-horvathskin_xb
summarize horvathskin_res

* residuals versus fitted values 
scatter horvathskin_res horvathskin_xb

* residuals versus chronological age 
scatter horvathskin_res page

histogram horvathskin_res, normal 

* GRIM AGE

regress dnamgrimage page
predict grimage_xb, xb
gen grimage_res = dnamgrimage-grimage_xb
summarize grimage_res

* residuals versus fitted values 
scatter grimage_res grimage_xb

* residuals versus chronological age 
scatter grimage_res page

histogram grimage_res, normal 

* MPOA


gen mpoa_dnamage = mpoa*page 
summarize mpoa_dnamage

* MPOA

regress mpoa_dnamage page
predict mpoa_xb, xb
gen mpoa_res = mpoa_dnamage-mpoa_xb
summarize mpoa_res

* residuals versus fitted values 
scatter mpoa_res mpoa_xb

* residuals versus chronological age 
scatter mpoa_res page

histogram mpoa_res, normal 




**************************************
**** INCLUSION VS EXCLUSION TABLE ****
**************************************



tabout raceth_f sex2 edubin_low dementia_lw  missing_2 using table_M_7-15-22.xls, ///
c(col) f(1) clab(Col_%) stats(chi2) nlab(Sample Size) ///
npos(row) percent ///
replace  


foreach var in page grimage_res mpoa_res  levine_res horthvath_res hannum_res horvathskin_res granulocytes monocytes lynphocytes {
asdoc ttest `var', by(missing_2) unequal rowappend stat(obs mean df t p) title(T-test) 
}

summarize page grimage_res mpoa_res  levine_res horthvath_res hannum_res horvathskin_res granulocytes monocytes lynphocytes



 






























