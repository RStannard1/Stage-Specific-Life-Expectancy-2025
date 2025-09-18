clear all

********************************************************************************
*** INSTALL PACKAGES ***
/*
ssc install rcsgen
ssc install gensplines // change from rcsgen to gensplines
ssc install stpm3
ssc install violinplot
ssc install stpp
ssc install standsurv
*/
********************************************************************************
*** SET WORKING DIRECTORY ***

// change this to your own local working directory
cd "R:/CanStat2/res53/Life_Expectancy/V2/Example code"

********************************************************************************
*** POPULATION MORTALITY FILE PREPARATION ***

use "https://pclambert.net/data/popmort.dta", clear
// alter to match data used for real analysis
replace _year=_year+26
save popmort.dta, replace

********************************************************************************
*** SET UP AND LOAD DATA ***

// imputations reduced from 50 to 3 for illustrative purposes
global nimp=3

use https://www.pclambert.net/data/colon, clear

drop subsite
drop year8594
drop agegrp

// change calendar year range in example dataset to match real data and popmort_ONS
replace yydx=yydx+26
replace exit=exit+26*365.25
replace dx=dx+26*365.25

// recode missing stage values = .
replace stage=. if stage==0
// stage in this example is Localised/Regional/Distant rather than stage I-IV
levelsof stage
global stageN=`r(r)'
di $stageN
tab stage, missing

// example dataset does not have deprivation information
// assign deprivation randomly and uniformly
gen dep_rand = runiform(0,100)
gen deprivation=.
replace deprivation=1 if dep_rand<20
replace deprivation=2 if dep_rand>=20 & dep_rand<40
replace deprivation=3 if dep_rand>=40 & dep_rand<60
replace deprivation=4 if dep_rand>=60 & dep_rand<80
replace deprivation=5 if dep_rand>=80

save colon.dta, replace

/*
// 6 month adjustment - not needed for example dataset but included for transparency
gen age_original=age // to keep an original copy
replace age=age+0.5
*/

*** WINSORIZE AGE AND CREATE SPLINE VARIABLES - FOR IMPUTATION ***
// ONE FOR MAIN EFFECTS AND ONE FOR INTERACTION TERMS/TIME-DEPENDENT EFFECTS
gensplines age, df(4) center(70) winsor(1 99) type(ns) gen(agercs)
gensplines age, df(2) center(70) winsor(1 99) type(ns) gen(age_rcs_int)

*** REFORMAT SEX VARIABLE ***
// labels already defined in this dataset
label define sexlbl 1 "Male" 2 "Female"
label values sex sexlbl

gen fem = sex==2 

********************************************************************************
*** STSET ***

stset exit                                                 					 ///
	,                                                        				 ///
	origin(dx)                                             					 ///
	failure(status==1,2)                                           			 ///
	scale(365.25)                                                			 ///
	id(id)                                                   			     ///
	enter(time mdy(1,1,2017))                               	 			 ///
	exit(time min(mdy(02,29,2020), dx+10.1*365.25))
	
keep if _st==1

su age if yydx>=2017 & yydx<=2019, detail
ret list
global age_med=`r(p50)'

su age if sex==1 & yydx>=2017 & yydx<=2019, detail
ret list
capture global age_med_male=`r(p50)'

su age if sex==2 & yydx>=2017 & yydx<=2019, detail
ret list
capture global age_med_fem=`r(p50)'

frame copy default age_hist, replace 

********************************************************************************
*** IMPUTE MISSING STAGE ***

sts gen H = na // Nelson-Aalen estimate of cumulative hazard

sts gen d = d // number of failures
sts gen n = n // number at risk

gen rev_d = 1-_d // one minus event indicator: 0 if dead, 1 if alive
// so that I can sort by deaths first

bysort _t (rev_d): gen first = _n==1 if _d==1 // gen first=1 for one person (who died) with each event time

gen td = _t*d/n if first & _d==1 // calculate Nelson-Aalen like estimator contribution for each event time

gen H1 = sum(td) // sum td to calculate Nelson-Aalen like estimator described by Keogh and Morris at each event time

// imputation interactions
gen lnt = log(_t)
rcsgen lnt, df(2) gen(lnt_rcs) if2(_d==1) // rcs for log(time)

forvalues i=1/2 {
	gen _dt`i' =  _d*lnt_rcs`i' // interaction between event indicator and log(time)
}

forvalues i=1/2 {
	gen H_age_rcs_int`i' = age_rcs_int`i'*H // interaction between age and Nelson-Aalen
}

forvalues i=1/2 {
	gen H1_age_rcs_int`i' = age_rcs_int`i'*H1 // interaction between age and H1
}

gen fem_H = fem*H // interaction between sex and Nelson-Aalen
gen fem_H1 = fem*H1 // interaction between sex and NA-like

/*
// copy yydx and change 2006 to 2007 for imputation only (very few obs in 2006 due to dx+10.1 in stset)
gen yydx_imp=yydx
replace yydx_imp=2007 if yydx==2006
*/

mi set flong

mi register imputed stage

gen fem_age1 = fem*age_rcs_int1
gen fem_age2 = fem*age_rcs_int2

set seed 18532

mi register regular fem agercs? yydx fem_age* deprivation H _d _dt* H1   	 ///
	H_age_rcs_int1 H_age_rcs_int2 H1_age_rcs_int1 H1_age_rcs_int2 fem_H fem_H1

mi impute mlogit stage = i.fem c.agercs? i.yydx fem_age*           			 ///
	i.deprivation i.yydx#(c.age_rcs_int? i.fem i.deprivation)			 	 ///
	H _d _dt* H1  											 				 ///
	H_age_rcs_int1 H_age_rcs_int2 H1_age_rcs_int1 H1_age_rcs_int2 	 		 ///
	fem_H fem_H1														   	 ///
		,                                                                  	 ///
		add(${nimp}) dots chaindots augment	
		
********************************************************************************
*** MERGE WITH POPMORT ***

gen _year = floor(yydx + _t)
gen _age = min(99, floor(age + _t))

sort sex _year _age
merge m:1 sex _year _age using "popmort.dta", keep(match master) keepusing(rate)
		
********************************************************************************
*** LOOK AT STAGE DISTRIBUTION IN EACH DATASET ***

save "imputed_data", replace	
		
********************************************************************************
*** FIT FPM TO EACH IMPUTED DATASET ***	
*** DEFINE CANDIDATE MODELS AND FIX KNOTS (# AND POSITION) ***

// this section fits models to original CC dataset (before imputation)
// here stage I = Localised, stage II = Regional and stage III = Distant

// loop for stage `s'

// winsor point 1%
// model 1

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " interaction model with age and sex tde"
	
	local s`s'_min = 1

	// interaction
	capture noisily stpm3 @ns(age, df(4) winsor(1 99)) ///
		@ns(age, df(2) winsor(1 99))#1.fem  i.fem      ///
		if stage==`s' & _mi_m==0                       ///
			,                                          ///
			bhazard(rate)                              ///
			scale(lncumhazard)                         ///
			df(5)                             		   ///
			tvc(@ns(age, df(2) winsor(1 99)) i.fem)    ///
			dftvc(2)                       			   ///
			initvaluesloop
			
	if _rc!=0 {
		local s`s'_min = 2
	}
	
	// baseline allknots()
	local base_knots1_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots1_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots1_s`s' = e(ef_age_knots1)
	// tde/int effects age allknots()
	local age_tde_knots1_s`s' = e(ef_age_knots2)	
	
	// winsor points for age spline
	local age_winsor1_`s' = e(ef_age_winsor1)
	// winsor points for tde/int age spline
	local age_tde_winsor1_`s' = e(ef_age_winsor2)
	
	local m1_`s' @ns(age, allknots(`age_knots1_s`s'') winsor(`age_winsor1_`s'', values)) @ns(age, allknots(`age_tde_knots1_s`s'') winsor(`age_tde_winsor1_`s'', values))#1.fem  i.fem
	local m1tvc_`s' tvc(@ns(age, allknots(`age_tde_knots1_s`s'') winsor(`age_tde_winsor1_`s'', values)) i.fem) allknotstvc(`tvc_knots1_s`s'', lntime)
	
}

// model 2

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " interaction model with age tde"

	// interaction
	capture noisily stpm3 @ns(age, df(4) winsor(1 99)) ///
		@ns(age, df(2) winsor(1 99))#1.fem  i.fem      ///
		if stage==`s' & _mi_m==0                       ///
			,                                          ///
			bhazard(rate)                              ///
			scale(lncumhazard)                         ///
			df(5)                             		   ///
			tvc(@ns(age, df(2) winsor(1 99)))          ///
			dftvc(2)                       			   ///
			initvaluesloop
			
	if _rc!=0 {
		local s`s'_min = 3
	}
	
	// baseline allknots()
	local base_knots2_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots2_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots2_s`s' = e(ef_age_knots1)
	// tde/int effects age allknots()
	local age_tde_knots2_s`s' = e(ef_age_knots2)	
	
	// winsor points for age spline
	local age_winsor2_`s' = e(ef_age_winsor1)
	// winsor points for tde/int age spline
	local age_tde_winsor2_`s' = e(ef_age_winsor2)
	
	local m2_`s' @ns(age, allknots(`age_knots2_s`s'') winsor(`age_winsor2_`s'', values)) @ns(age, allknots(`age_tde_knots2_s`s'') winsor(`age_tde_winsor2_`s'', values))#1.fem  i.fem
	local m2tvc_`s' tvc(@ns(age, allknots(`age_tde_knots2_s`s'') winsor(`age_tde_winsor2_`s'', values))) allknotstvc(`tvc_knots2_s`s'', lntime)
	
}

// 2% winsor

// model 3

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " interaction model"

	// interaction
	capture noisily stpm3 @ns(age, df(4) winsor(2 98)) ///
		@ns(age, df(2) winsor(2 98))#1.fem  i.fem      ///
		if stage==`s' & _mi_m==0                       ///
			,                                          ///
			bhazard(rate)                              ///
			scale(lncumhazard)                         ///
			df(5)                             		   ///
			tvc(@ns(age, df(2) winsor(2 98)))          ///
			dftvc(2)                       			   ///
			initvaluesloop
			
	if _rc!=0 {
		local s`s'_min = 4
	}
	
	// baseline allknots()
	local base_knots3_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots3_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots3_s`s' = e(ef_age_knots1)
	// tde/int effects age allknots()
	local age_tde_knots3_s`s' = e(ef_age_knots2)	
	
	// winsor points for age spline
	local age_winsor3_`s' = e(ef_age_winsor1)
	// winsor points for tde/int age spline
	local age_tde_winsor3_`s' = e(ef_age_winsor2)
	
	local m3_`s' @ns(age, allknots(`age_knots3_s`s'') winsor(`age_winsor3_`s'', values)) @ns(age, allknots(`age_tde_knots3_s`s'') winsor(`age_tde_winsor3_`s'', values))#1.fem  i.fem
	local m3tvc_`s' tvc(@ns(age, allknots(`age_tde_knots3_s`s'') winsor(`age_tde_winsor3_`s'', values))) allknotstvc(`tvc_knots3_s`s'', lntime)
	
}

// model 4

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " no interaction model"
	
	capture noisily stpm3 @ns(age, df(4) winsor(2 98)) i.fem ///
		if stage==`s' & _mi_m==0                             ///
		,                                                    ///
		bhazard(rate)                                        ///
		scale(lncumhazard)                                   ///
		df(5)                                       		 ///
		tvc(@ns(age, df(2) winsor(2 98)))                    ///
		dftvc(2)                                 			 ///
		initvaluesloop		
		
	if _rc!=0 {
		local s`s'_min = 5 
	}
	
	// baseline allknots()
	local base_knots4_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots4_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots4_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots4_s`s' = e(ef_age_knots2)	
	
	// winsor points for age spline
	local age_winsor4_`s' = e(ef_age_winsor1)
	// winsor points for tde age spline
	local age_tde_winsor4_`s' = e(ef_age_winsor2)
	
	local m4_`s' @ns(age, allknots(`age_knots4_s`s'') winsor(`age_winsor4_`s'', values)) i.fem
	local m4tvc_`s' tvc(@ns(age, allknots(`age_tde_knots4_s`s'') winsor(`age_tde_winsor4_`s'', values))) allknotstvc(`tvc_knots4_s`s'', lntime)
	
}

// model 5

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " no time-dependent effects model"
	
	capture noisily stpm3 @ns(age, df(4) winsor(2 98)) ///
		@ns(age, df(2) winsor(2 98))#1.fem i.fem       ///
		if stage==`s' & _mi_m==0                       ///
			,                                          ///
			bhazard(rate)                              ///
			scale(lncumhazard)                         ///
			df(5)                             		   ///
			initvaluesloop	
			
	if _rc!=0 {
		local s`s'_min = 6
	}		
	
	// baseline allknots()
	local base_knots5_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots5_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots5_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots5_s`s' = e(ef_age_knots2)
	
	// winsor points for age spline
	local age_winsor5_`s' = e(ef_age_winsor1)
	
	local m5_`s' @ns(age, allknots(`age_knots5_s`s'') winsor(`age_winsor5_`s'', values)) @ns(age, allknots(`age_tde_knots5_s`s'') winsor(`age_winsor5_`s'', values))#1.fem  i.fem
	local m5tvc_`s' ""
	
}

// model 6

forvalues s=1/$stageN {	
	
	di "stage " "`s'" " no time-dependent effects or interactions model"
	
	capture noisily stpm3 @ns(age, df(4) winsor(2 98)) i.fem ///
		if stage==`s' & _mi_m==0                             ///
			,                                                ///
			bhazard(rate)                                    ///
			scale(lncumhazard)                               ///
			df(5)                                   		 ///
			initvaluesloop	
			
	if _rc!=0 {
		local s`s'_min = 7
	}
	
	// baseline allknots()
	local base_knots6_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots6_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots6_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots6_s`s' = e(ef_age_knots2)
	
	// winsor points for age spline
	local age_winsor6_`s' = e(ef_age_winsor1)
	
	local m6_`s' @ns(age, allknots(`age_knots6_s`s'') winsor(`age_winsor6_`s'', values)) i.fem
	local m6tvc_`s' ""
}

// 4 baseline df (instead of 5)

// model 7

forvalues s=1/$stageN {	
	
	di "less baseline df (4)"
	
	capture noisily stpm3 @ns(age, df(4) winsor(2 98)) i.fem ///
		if stage==`s' & _mi_m==0                             ///
			,                                                ///
			bhazard(rate)                                    ///
			scale(lncumhazard)                               ///
			df(4)                                  			 ///
			initvaluesloop		  
			
	if _rc!=0 {
		local s`s'_min = 8
	}
	
	// baseline allknots()
	local base_knots7_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots7_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots7_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots7_s`s' = e(ef_age_knots2)
	
	// winsor points for age spline
	local age_winsor7_`s' = e(ef_age_winsor1)
	
	local m7_`s' @ns(age, allknots(`age_knots7_s`s'') winsor(`age_winsor7_`s'', values)) i.fem
	local m7tvc_`s' ""
}

// model 8

forvalues s=1/$stageN {	
	
	di "less age df"
	
	capture noisily stpm3 @ns(age, df(3) winsor(2 98)) i.fem ///
		if stage==`s' & _mi_m==0                             ///
			,                                                ///
			bhazard(rate)                                    ///
			scale(lncumhazard)                               ///
			df(4)                                  			 ///
			initvaluesloop		                                                
	
	if _rc!=0 {
		local s`s'_min = 9
	}
	
	// baseline allknots()
	local base_knots8_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots8_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots8_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots8_s`s' = e(ef_age_knots2)
	
	// winsor points for age spline
	local age_winsor8_`s' = e(ef_age_winsor1)
	
	local m8_`s' @ns(age, allknots(`age_knots8_s`s'') winsor(`age_winsor8_`s'', values)) i.fem
	local m8tvc_`s' ""
}

// model 9
// this model is added purely for s1 in this example dataset

forvalues s=1/$stageN {	
	
	di "less baseline df (3)"
	
	capture noisily stpm3 @ns(age, df(3) winsor(2 98)) i.fem ///
		if stage==`s' & _mi_m==0                             ///
			,                                                ///
			bhazard(rate)                                    ///
			scale(lncumhazard)                               ///
			df(3)                                  			 ///
			initvaluesloop		                                                
	
	// baseline allknots()
	local base_knots9_s`s' = e(knots)
	// tde allknotstvc()
	local tvc_knots9_s`s' = e(knots_tvc)	
	// main effects age allknots()
	local age_knots9_s`s' = e(ef_age_knots1)
	// tde effects age allknots()
	local age_tde_knots9_s`s' = e(ef_age_knots2)
	
	// winsor points for age spline
	local age_winsor9_`s' = e(ef_age_winsor1)
	
	local m9_`s' @ns(age, allknots(`age_knots9_s`s'') winsor(`age_winsor9_`s'', values)) i.fem
	local m9tvc_`s' ""
}
	 
*** FIT MODELS - for each impuated dataset and each stage ***
forvalues s=1/$stageN {
	forvalues j=1/$nimp {
		cap erase stpm3_s`s'_m`j'.ster
	}
}

cap frame drop model_check
frame create model_check imp_number stage model_number str8 winsor_point

forvalues s=1/$stageN {
	di `s`s'_min'
}

forvalues j=1/$nimp {
	
	preserve 
	
	mi extract `j', clear
	
	forvalues s=1/$stageN {
		
		di "imputation " "`j'" " stage " "`s'"
		
		forvalues i=`s`s'_min'/9 {
			
			di "model " `i'
			
			capture noisily stpm3 `m`i'_`s'' if stage==`s' 					 ///
				,                                          					 ///
				bhazard(rate)                              					 ///
				scale(lncumhazard)                         					 ///
				allknots(`base_knots`i'_s`s'', lntime)     					 ///
				`m`i'tvc_`s''                              					 ///
				initvaluesloop
				
			if "`e(converged)'"=="1" {
				frame post model_check (`j') (`s') (`i') ("`age_winsor`i'_`s''")
				
				est store stpm3_s`s'_m`j'
				est save stpm3_s`s'_m`j', replace
				
				continue, break
				
			}
			
		}
		
	}	
	
	restore
}	

frame model_check: save "model_check", replace
	
********************************************************************************
*** ESTIMATE LIFE EXPECTANCY ***

cap frame drop life_exp
frame create life_exp
frame change life_exp

range age 30 90 61
gen fem=.
gen sex=.
gen tt=90

forvalues k=1/$nimp {
	
	forvalues s=1/$stageN {
		
		est use stpm3_s`s'_m`k'
	
		cap drop es_m es_f
	
		predict le_m_s`s'_`k' le_f_s`s'_`k'                         		 ///
			,                                                       		 ///
			rmst ci merge                                           		 ///
			at1(fem 0, obsvalues)                                   		 ///
			at2(fem 1, obsvalues)                                   		 ///
			timevar(tt)                                             		 ///
			expsurv(using("popmort") 		 								 ///
				datediag(2019-1-1)                                  		 ///
				agediag(age)                                        		 ///
				pmrate(rate)			                            		 ///
				pmage(_age)				                            		 ///
				pmyear(_year)                                       		 ///
				pmother(sex)                                        		 ///
				at1(sex 1)                                          		 ///
				at2(sex 2)                                          		 ///
				pmmaxyear(2019)                                     		 ///
				expvars(es_m es_f))
	}
}

save "life_exp_mi", replace

// take the log of the le estimates and appoximate within-imputation variance
forvalues s=1/$stageN {
	forvalues i=1/$nimp {
		gen ln_le_m_s`s'_`i' = ln(le_m_s`s'_`i')
		gen ln_le_f_s`s'_`i' = ln(le_f_s`s'_`i')
		
		// Approximate variance on log scale.
		gen var_ln_le_m_s`s'_`i' = ///
			((ln(le_m_s`s'_`i'_uci)-ln(le_m_s`s'_`i'))/1.96)^2
		gen var_ln_le_f_s`s'_`i' = ///
			((ln(le_f_s`s'_`i'_uci)-ln(le_f_s`s'_`i'))/1.96)^2
	}
}

drop le_*

// combine using Rubin's Rules	
forvalues s=1/$stageN {
	// take the mean and exponentiate
	egen ln_le_m_s`s' = rowmean(ln_le_m_s`s'*)
	gen le_m_s`s' = exp(ln_le_m_s`s')
	
	egen ln_le_f_s`s' = rowmean(ln_le_f_s`s'*)
	gen le_f_s`s' = exp(ln_le_f_s`s')
	
	// calculate se and CIs
	egen ln_le_m_s`s'_B = rowsd(ln_le_m_s`s'_*) 
	egen ln_le_m_s`s'_U = rowmean(var_ln_le_m_s`s'_*)
	gen ln_le_m_s`s'_se = sqrt(ln_le_m_s`s'_U+(1+1/${nimp})*ln_le_m_s`s'_B^2)
	gen le_m_s`s'_lci = exp(ln_le_m_s`s'-1.96*ln_le_m_s`s'_se)
	gen le_m_s`s'_uci = exp(ln_le_m_s`s'+1.96*ln_le_m_s`s'_se)
	
	egen ln_le_f_s`s'_B = rowsd(ln_le_f_s`s'_*) 
	egen ln_le_f_s`s'_U = rowmean(var_ln_le_f_s`s'_*)
	gen ln_le_f_s`s'_se = sqrt(ln_le_f_s`s'_U+(1+1/${nimp})*ln_le_f_s`s'_B^2)
	gen le_f_s`s'_lci = exp(ln_le_f_s`s'-1.96*ln_le_f_s`s'_se)
	gen le_f_s`s'_uci = exp(ln_le_f_s`s'+1.96*ln_le_f_s`s'_se)
}

// generate attained life expectancy
foreach sex in "m" "f" {
	forvalues s=1/$stageN {
		gen lifeexp_`sex'_s`s' = le_`sex'_s`s'+age
		
		gen lifeexp_`sex'_s`s'_lci = le_`sex'_s`s'_lci+age
		gen lifeexp_`sex'_s`s'_uci = le_`sex'_s`s'_uci+age
	}
}

// generate expected attained life expectancy (for comparison)
foreach sex in "m" "f" {
	gen pop_le_`sex' =es_`sex'+age
}

save "life_exp", replace

********************************************************************************
*** PLOTS ***

// overall median
local row_n=${age_med}-29

global pop_le_m_med=pop_le_m[`row_n']
global s1_le_m_med=lifeexp_m_s1[`row_n']
global s2_le_m_med=lifeexp_m_s2[`row_n']
global s3_le_m_med=lifeexp_m_s3[`row_n']

global pop_le_f_med=pop_le_f[`row_n']
global s1_le_f_med=lifeexp_f_s1[`row_n']
global s2_le_f_med=lifeexp_f_s2[`row_n']
global s3_le_f_med=lifeexp_f_s3[`row_n']

global pop_le_m_med: di %3.1f ${pop_le_m_med}
global s1_le_m_med: di %3.1f ${s1_le_m_med}
global s2_le_m_med: di %3.1f ${s2_le_m_med}
global s3_le_m_med: di %3.1f ${s3_le_m_med}

global pop_le_f_med: di %3.1f ${pop_le_f_med}
global s1_le_f_med: di %3.1f ${s1_le_f_med}
global s2_le_f_med: di %3.1f ${s2_le_f_med}
global s3_le_f_med: di %3.1f ${s3_le_f_med}

// sex-specific medians
local row_n_m=${age_med_male}-29
global pop_le_m_med2=pop_le_m[`row_n_m']
global s1_le_m_med2=lifeexp_m_s1[`row_n_m']
global s2_le_m_med2=lifeexp_m_s2[`row_n_m']
global s3_le_m_med2=lifeexp_m_s3[`row_n_m']

local row_n_f=${age_med_fem}-29
global pop_le_f_med2=pop_le_f[`row_n_f']
global s1_le_f_med2=lifeexp_f_s1[`row_n_f']
global s2_le_f_med2=lifeexp_f_s2[`row_n_f']
global s3_le_f_med2=lifeexp_f_s3[`row_n_f']

global pop_le_m_med2: di %3.1f ${pop_le_m_med2}
global s1_le_m_med2: di %3.1f ${s1_le_m_med2}
global s2_le_m_med2: di %3.1f ${s2_le_m_med2}
global s3_le_m_med2: di %3.1f ${s3_le_m_med2}

global pop_le_f_med2: di %3.1f ${pop_le_f_med2}
global s1_le_f_med2: di %3.1f ${s1_le_f_med2}
global s2_le_f_med2: di %3.1f ${s2_le_f_med2}
global s3_le_f_med2: di %3.1f ${s3_le_f_med2}

*** FOR PANEL PLOTS ***  

// sex panels
twoway (line pop_le_m age 													 ///
			,                                                                ///
			lpattern(dash)                            						 ///
			lcolor(black))            										 ///
		(line lifeexp_m_s1 lifeexp_m_s2 lifeexp_m_s3 age 		 			 ///
			if age>=40 														 ///
			,                                                                ///
			lpattern(solid solid solid)                            	 		 ///
			lcolor("173 40 113" lavender "49 167 173"))            	 		 ///
		(area age age                                                        ///
			,                                                                ///
			sort                                                             ///
			color(gray)                                                      ///
			fintensity(30))                                                  ///
		(function y=${age_med_male}, horizontal range(30 100) lcolor(gray))  ///
			,                                                                ///
			ytitle("Predicted Age at Death (years)", color(none))            ///
			xtitle("")                   									 ///
			legend(off)                                                      ///
			fxsize(100)                        			  					 ///
			name(m, replace)                                        		 ///
			plotregion(margin(zero))	                                     ///
			text(58 89.5 "{bf:LE at age ${age_med_male}:}",        			 ///
				place(w) color(gray) size(small))                            ///
			text(53 89.5 "Population - ${pop_le_m_med2}",               	 ///
				place(w) color(black) size(small))                           ///
			text(48 89.5 "Stage I - ${s1_le_m_med2}",                        ///
				place(w) color("173 40 113") size(small))                    ///
			text(43 89.5 "Stage II - ${s2_le_m_med2}",                       ///
				place(w) color(lavender) size(small))                        ///
			text(38 89.5 "Stage III - ${s3_le_m_med2}",                      ///
				place(w) color("49 167 173") size(small))                    ///
			xscale(range(30 90))                                             ///
			xlabel(30(10)90)                                                 ///
			yscale(range(30 100))                                            ///
			ylabel(30(10)100, angle(h) format(%2.0f))                  
	
graph export "le_m.pdf", replace
graph export "le_m.png", replace
		
twoway (line pop_le_f age 													 ///
			,                                                                ///
			lpattern(dash)                            						 ///
			lcolor(black))            										 ///
		(line lifeexp_f_s1 lifeexp_f_s2 lifeexp_f_s3 age 		 			 ///
			if age>=40 														 ///
			,                                                                ///
			lpattern(solid solid solid)                            			 ///
			lcolor("173 40 113" lavender "49 167 173"))            	 		 ///
		(area age age                                                        ///
			,                                                                ///
			sort                                                             ///
			color(gray)                                                      ///
			fintensity(30))                                                  ///
		(function y=${age_med_fem}, horizontal range(30 100) lcolor(gray))   ///
			,                                                                ///
			ytitle("Predicted Age at Death (years)", color(none))            ///
			xtitle("")                   									 ///
			legend(off)                                                      ///
			fxsize(100)                        			  					 ///
			name(f, replace)                                       			 ///
			plotregion(margin(zero))	                                     ///
			text(58 89.5 "{bf:LE at age ${age_med_fem}:}",         			 ///
				place(w) color(gray) size(small))                            ///
			text(53 89.5 "Population - ${pop_le_f_med2}",               	 ///
				place(w) color(black) size(small))                           ///
			text(48 89.5 "Stage I - ${s1_le_f_med2}",                        ///
				place(w) color("173 40 113") size(small))                    ///
			text(43 89.5 "Stage II - ${s2_le_f_med2}",                       ///
				place(w) color(lavender) size(small))                        ///
			text(38 89.5 "Stage III - ${s3_le_f_med2}",                      ///
				place(w) color("49 167 173") size(small))                    ///
			xscale(range(30 90))                                             ///
			xlabel(30(10)90)                                                 ///
			yscale(range(30 100))                                            ///
			ylabel(30(10)100,angle(h) format(%2.0f))
				
graph export "le_f.pdf", replace
graph export "le_f.png", replace

use "imputed_data", clear

// male
preserve
	
keep if yydx>=2017 & yydx<=2019 & _mi_m!=0 & sex==1
	
gen one=1
collapse (sum) one, by(age stage)
	
replace one=floor(one/${nimp})
	
violinplot age [fw=one] 										 		 	 ///
	, 		             		  	   			 		 				 	 ///
	over(stage)									 		 				 	 ///
	horizontal									 		 				 	 ///
	overlay 									 		 				 	 ///
	left 										 		 				 	 ///
	pstyles(1 2 3)                             		 				 	 	 ///
	p1(color("173 40 113"))           			 		 				 	 ///
	p2(color(lavender))           				 		 				 	 ///
	p3(color("49 167 173"))           			 		 				 	 ///
	xtitle("")                         			 		 				 	 ///
	xscale(range(30 90) noline) 	   			 		 				 	 ///
	xlabel(30(10)90, labcolor(none) tlcolor(none) nogrid)				 	 ///
	range(40 90) 										 				 	 ///
	absolute											 				 	 ///
	dscale(.)											 				 	 ///
	yscale(range(0 30) noline)                	 				 	 		 ///
	ylabel(0 30, labcolor(none) tlcolor(none)nogrid format(%9.0f))           ///
	text(0 60 "Stage-specific age distribution", 						 	 ///
		size(small) placement(south)) 									 	 ///
	ytitle("this is the title", color(none))     		 				 	 ///
	name(violin_m, replace)      		 	 				 		 		 ///
	plotregion(margin(zero))           			 		 				 	 ///
	fxsize(100)                        			 		 				 	 ///
	fysize(15)                                   		 				 	 ///
	line(lwidth(vthin)) 								 				 	 ///
	nobox nowhiskers nomedian                            				 	 ///
	text(30 30 "Colon Example Data - Male", 						 	 	 ///
		size(medsmall) placement(se)) 									 	 ///
	legend(off)
		
graph combine violin_m m, col(1) imargin(0 0 0 0) saving("le_dens_m", replace)
	
		
graph export "le_dens_m.pdf", replace
graph export "le_dens_m.png", replace
	
restore

// female
preserve
	
keep if yydx>=2017 & yydx<=2019 & _mi_m!=0 & sex==2
	
gen one=1
collapse (sum) one, by(age stage)
	
replace one=floor(one/${nimp})
	
violinplot age [fw=one] 										 		 	 ///
	, 		             		  	   			 		 				 	 ///
	over(stage)									 		 				 	 ///
	horizontal									 		 				 	 ///
	overlay 									 		 				 	 ///
	left 										 		 				 	 ///
	pstyles(1 2 3)                             		 				 	 	 ///
	p1(color("173 40 113"))           			 		 				 	 ///
	p2(color(lavender))           				 		 				 	 ///
	p3(color("49 167 173"))           			 		 				 	 ///
	xtitle("")                         			 		 				 	 ///
	xscale(range(30 90) noline) 	   			 		 				 	 ///
	xlabel(30(10)90, labcolor(none) tlcolor(none) nogrid)				 	 ///
	range(40 90) 										 				 	 ///
	absolute											 				 	 ///
	dscale(.)											 				 	 ///
	yscale(range(0 30) noline)                	 				 	 		 ///
	ylabel(0 30, labcolor(none) tlcolor(none)nogrid format(%9.0f))           ///
	text(0 60 "Stage-specific age distribution", 						 	 ///
		size(small) placement(south)) 									 	 ///
	ytitle("this is the title", color(none))     		 				 	 ///
	name(violin_f, replace)      		 	 				 		 		 ///
	plotregion(margin(zero))           			 		 				 	 ///
	fxsize(100)                        			 		 				 	 ///
	fysize(15)                                   		 				 	 ///
	line(lwidth(vthin)) 								 				 	 ///
	nobox nowhiskers nomedian                            				 	 ///
	text(30 30 "Colon Example Data - Female", 						 	 	 ///
		size(medsmall) placement(se)) 									 	 ///
	legend(off)
		
graph combine violin_f f, col(1) imargin(0 0 0 0) saving("le_dens_f", replace)
	
		
graph export "le_dens_f.pdf", replace
graph export "le_dens_f.png", replace
	
restore

********************************************************************************
*** COMPARE TO PP FOR VALIDATION ***

frame change default

use "imputed_data", clear

cap drop rs_pp

forvalues k=1/$nimp {
	cap frame drop rs_pp_`k'
}

recode age (min/45=1) (45/55=2) (55/65=3) (65/75=4) (75/max=5), gen(ICSSagegrp)

*** POHAR-PERME FOR VALIDATION ***

// estimate PP and reshape to desired wide format
forvalues k=1/$nimp {
	
	preserve
	
	mi extract `k', clear
	
	stpp rs_`k' using "popmort" 				 							 ///
		,                                                      				 ///
		agediag(age)                                           				 ///
		datediag(dx)                                     					 ///
		pmother(sex)                                           				 ///
		by(sex ICSSagegrp stage)                               				 ///
		list(1(1)10)                                           				 ///
		frame(rs_pp_`k')
	
	frame rs_pp_`k' {
		
		drop Natrisk
		
		rename PP rs_`k'
		rename PP_lci rs_`k'_lci
		rename PP_uci rs_`k'_uci
	
		reshape wide rs_`k' rs_`k'_lci rs_`k'_uci, 							 ///
			i(time ICSSagegrp stage) j(sex)
		
		rename rs_`k'1 rs_m_`k'
		rename rs_`k'2 rs_f_`k'
		rename rs_`k'_lci1 rs_m_`k'_lci
		rename rs_`k'_lci2 rs_f_`k'_lci
		rename rs_`k'_uci1 rs_m_`k'_uci
		rename rs_`k'_uci2 rs_f_`k'_uci
		
		reshape wide rs_m_`k' rs_f_`k' rs_m_`k'_lci rs_f_`k'_lci 			 ///
			rs_m_`k'_uci rs_f_`k'_uci, 							 			 ///
				i(time stage) j(ICSSagegrp)
	}
	
	foreach sex in m f {
		forvalues i=1/5 {
			frame rs_pp_`k' {
				
				rename rs_`sex'_`k'`i' rs_`sex'_a`i'_`k'
				rename rs_`sex'_`k'_lci`i' rs_`sex'_a`i'_`k'_lci
				rename rs_`sex'_`k'_uci`i' rs_`sex'_a`i'_`k'_uci
				
			}
		}
	}
	
	foreach sex in m f {
		forvalues i=1/5 {
			
			local varnames`k' `varnames`k'' rs_`sex'_a`i'_`k'
			local varnames`k'_lci `varnames`k'_lci' rs_`sex'_a`i'_`k'_lci
			local varnames`k'_uci `varnames`k'_uci' rs_`sex'_a`i'_`k'_uci
			
		}
	}
	
	frame rs_pp_`k': reshape wide 				   				     		 ///
		`varnames`k'' `varnames`k'_lci' `varnames`k'_uci', 				     ///
			i(time) j(stage)
	
	foreach sex in m f {
		forvalues i=1/5 {
			forvalues s=1/$stageN {
				frame rs_pp_`k' {
					
					rename rs_`sex'_a`i'_`k'`s' rs_`sex'_a`i'_s`s'_`k'
					rename rs_`sex'_a`i'_`k'_lci`s' rs_`sex'_a`i'_s`s'_`k'_lci
					rename rs_`sex'_a`i'_`k'_uci`s' rs_`sex'_a`i'_s`s'_`k'_uci
					
				}
			}
		}
	}
	
	restore
	
}

// merge data frames
cap frame drop rs_pp
frame copy rs_pp_1 rs_pp_temp, replace
frame change rs_pp_temp

forvalues k=2/$nimp {
	frlink m:1 time, frame(rs_pp_`k')
}

// age group 1
forvalues k=2/$nimp {	
	
	foreach sex in m f {

		forvalues s=1/$stageN {	
				
			local getnames_`k'_a1 `getnames_`k'_a1' rs_`sex'_a1_s`s'_`k'
			local getnames_`k'_a1_lci `getnames_`k'_a1_lci' rs_`sex'_a1_s`s'_`k'_lci
			local getnames_`k'_a1_uci `getnames_`k'_a1_uci' rs_`sex'_a1_s`s'_`k'_uci
		
		}
		
	}
	
	frget `getnames_`k'_a1' `getnames_`k'_a1_lci' `getnames_`k'_a1_uci', 	 ///
		from(rs_pp_`k')
	
}		

preserve
keep rs_m_a1_* rs_f_a1_* time
frame copy rs_pp_temp rs_pp_a1, replace
restore

forvalues k=2/$nimp {
	drop `getnames_`k'_a1' `getnames_`k'_a1_lci' `getnames_`k'_a1_uci'
}

// age group 2
forvalues k=2/$nimp {	
	
	foreach sex in m f {

		forvalues s=1/$stageN {	
				
			local getnames_`k'_a2 `getnames_`k'_a2' rs_`sex'_a2_s`s'_`k'
			local getnames_`k'_a2_lci `getnames_`k'_a2_lci' rs_`sex'_a2_s`s'_`k'_lci
			local getnames_`k'_a2_uci `getnames_`k'_a2_uci' rs_`sex'_a2_s`s'_`k'_uci
		
		}
		
	}
	
	frget `getnames_`k'_a2' `getnames_`k'_a2_lci' `getnames_`k'_a2_uci',  	 ///
		from(rs_pp_`k')
	
}		

preserve
keep rs_m_a2_* rs_f_a2_* time
frame copy rs_pp_temp rs_pp_a2, replace
restore

forvalues k=2/$nimp {
	drop `getnames_`k'_a2' `getnames_`k'_a2_lci' `getnames_`k'_a2_uci'
}

// age group 3
forvalues k=2/$nimp {	
	
	foreach sex in m f {

		forvalues s=1/$stageN {	
				
			local getnames_`k'_a3 `getnames_`k'_a3' rs_`sex'_a3_s`s'_`k'
			local getnames_`k'_a3_lci `getnames_`k'_a3_lci' rs_`sex'_a3_s`s'_`k'_lci
			local getnames_`k'_a3_uci `getnames_`k'_a3_uci' rs_`sex'_a3_s`s'_`k'_uci
		
		}
		
	}
	
	frget `getnames_`k'_a3' `getnames_`k'_a3_lci' `getnames_`k'_a3_uci', 	 ///
		from(rs_pp_`k')
	
}		

preserve
keep rs_m_a3_* rs_f_a3_* time
frame copy rs_pp_temp rs_pp_a3, replace
restore

forvalues k=2/$nimp {
	drop `getnames_`k'_a3' `getnames_`k'_a3_lci' `getnames_`k'_a3_uci'
}

// age group 4
forvalues k=2/$nimp {	
	
	foreach sex in m f {

		forvalues s=1/$stageN {	
				
			local getnames_`k'_a4 `getnames_`k'_a4' rs_`sex'_a4_s`s'_`k'
			local getnames_`k'_a4_lci `getnames_`k'_a4_lci' rs_`sex'_a4_s`s'_`k'_lci
			local getnames_`k'_a4_uci `getnames_`k'_a4_uci' rs_`sex'_a4_s`s'_`k'_uci
		
		}
		
	}
	
	frget `getnames_`k'_a4' `getnames_`k'_a4_lci' `getnames_`k'_a4_uci', 	 ///
		from(rs_pp_`k')
	
}		

preserve
keep rs_m_a4_* rs_f_a4_* time
frame copy rs_pp_temp rs_pp_a4, replace
restore

forvalues k=2/$nimp {
	drop `getnames_`k'_a4' `getnames_`k'_a4_lci' `getnames_`k'_a4_uci'
}

// age group 5
forvalues k=2/$nimp {	
	
	foreach sex in m f {

		forvalues s=1/$stageN {	
				
			local getnames_`k'_a5 `getnames_`k'_a5' rs_`sex'_a5_s`s'_`k'
			local getnames_`k'_a5_lci `getnames_`k'_a5_lci' rs_`sex'_a5_s`s'_`k'_lci
			local getnames_`k'_a5_uci `getnames_`k'_a5_uci' rs_`sex'_a5_s`s'_`k'_uci
		
		}
		
	}
	
	frget `getnames_`k'_a5' `getnames_`k'_a5_lci' `getnames_`k'_a5_uci', 	 ///
		from(rs_pp_`k')
	
}		

preserve
keep rs_m_a5_* rs_f_a5_* time
frame copy rs_pp_temp rs_pp_a5, replace
restore

forvalues k=2/$nimp {
	drop `getnames_`k'_a5' `getnames_`k'_a5_lci' `getnames_`k'_a5_uci'
}

// rubin's rules
// males first (to not exceed maximum var number)
forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		forvalues k=1/$nimp {
			
			gen ln_rs_m_a`i'_s`s'_`k' = ln(rs_m_a`i'_s`s'_`k')
			
			// Approximate variance on log scale
			gen var_ln_rs_m_a`i'_s`s'_`k' = ///
				((ln(rs_m_a`i'_s`s'_`k'_uci)-ln(rs_m_a`i'_s`s'_`k'))/1.96)^2
		}
	}
}

forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		
		// take the mean and exponentiate
		egen ln_rs_m_a`i'_s`s' = rowmean(ln_rs_m_a`i'_s`s'*)
		gen rs_m_a`i'_s`s' = exp(ln_rs_m_a`i'_s`s')
		
		// calculate se and CIs
		egen ln_rs_m_a`i'_s`s'_B = rowsd(ln_rs_m_a`i'_s`s'_*) 
		egen ln_rs_m_a`i'_s`s'_U = rowmean(var_ln_rs_m_a`i'_s`s'_*)
		gen ln_rs_m_a`i'_s`s'_se = ///
			sqrt(ln_rs_m_a`i'_s`s'_U+(1+1/${nimp})*ln_rs_m_a`i'_s`s'_B^2)
		gen rs_m_a`i'_s`s'_lci = exp(ln_rs_m_a`i'_s`s'-1.96*ln_rs_m_a`i'_s`s'_se)
		gen rs_m_a`i'_s`s'_uci = exp(ln_rs_m_a`i'_s`s'+1.96*ln_rs_m_a`i'_s`s'_se)
	
	}
}

forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		
		drop ln_rs_m_a`i'_s`s' ln_rs_m_a`i'_s`s'_B ln_rs_m_a`i'_s`s'_U ln_rs_m_a`i'_s`s'_se
		
		forvalues k=1/$nimp {
		
			drop rs_m_a`i'_s`s'_`k' rs_m_a`i'_s`s'_`k'_uci rs_m_a`i'_s`s'_`k'_lci ln_rs_m_a`i'_s`s'_`k' var_ln_rs_m_a`i'_s`s'_`k'
			
		}
	}
}

// females
forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		forvalues k=1/$nimp {
			
			gen ln_rs_f_a`i'_s`s'_`k' = ln(rs_f_a`i'_s`s'_`k')
			
			// Approximate variance on log scale
			gen var_ln_rs_f_a`i'_s`s'_`k' = ///
				((ln(rs_f_a`i'_s`s'_`k'_uci)-ln(rs_f_a`i'_s`s'_`k'))/1.96)^2
		}
	}
}

forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		
		// take the mean and exponentiate
		egen ln_rs_f_a`i'_s`s' = rowmean(ln_rs_f_a`i'_s`s'*)
		gen rs_f_a`i'_s`s' = exp(ln_rs_f_a`i'_s`s')
		
		// calculate se and CIs
		egen ln_rs_f_a`i'_s`s'_B = rowsd(ln_rs_f_a`i'_s`s'_*) 
		egen ln_rs_f_a`i'_s`s'_U = rowmean(var_ln_rs_f_a`i'_s`s'_*)
		gen ln_rs_f_a`i'_s`s'_se = ///
			sqrt(ln_rs_f_a`i'_s`s'_U+(1+1/${nimp})*ln_rs_f_a`i'_s`s'_B^2)
		gen rs_f_a`i'_s`s'_lci = exp(ln_rs_f_a`i'_s`s'-1.96*ln_rs_f_a`i'_s`s'_se)
		gen rs_f_a`i'_s`s'_uci = exp(ln_rs_f_a`i'_s`s'+1.96*ln_rs_f_a`i'_s`s'_se)
	}
}

forvalues i=1/5 {
	
	frame change rs_pp_a`i'
	
	forvalues s=1/$stageN {
		
		drop ln_rs_f_a`i'_s`s' ln_rs_f_a`i'_s`s'_B ln_rs_f_a`i'_s`s'_U ln_rs_f_a`i'_s`s'_se
		
		forvalues k=1/$nimp {
		
			drop rs_f_a`i'_s`s'_`k' rs_f_a`i'_s`s'_`k'_uci rs_f_a`i'_s`s'_`k'_lci ln_rs_f_a`i'_s`s'_`k' var_ln_rs_f_a`i'_s`s'_`k'
			
		}
	}
}

// merge back into one frame
frame copy rs_pp_a1 rs_pp, replace
frame change rs_pp

forvalues i=2/5 {
	
	frlink m:1 time, frame(rs_pp_a`i')
	
}

forvalues i=2/5 {	
	
	frget *, from(rs_pp_a`i') exclude(time)
	
}

drop rs_pp_a*

save "rs_pp", replace

frame change default

*** MARGINAL RS ***

forvalues k=1/$nimp {
	forvalues s=1/$stageN {
		preserve
		
		keep if stage==`s' & _t0==0
		        
		est use stpm3_s`s'_m`k'
		
		range timevar 0 10 101
		
		foreach sex in m f {
			foreach agegrp in 1 2 3 4 5 {
				local atvars`s'`k' `atvars`s'`k'' rs_fpm_`sex'_a`agegrp'_s`s'_`k'
			}
		}
		
		di "s " `s' " 	k " `k'
		
		standsurv                  											 ///
			,                      											 ///
			ci					   											 ///
			surv                   											 ///
			atvars(`atvars`s'`k'') 											 ///
			over(fem ICSSagegrp)   											 ///
			timevar(timevar)       											 ///
			frame(rs_fpm_`s', mergecreate) 
		
		restore
	
	}
}

// split into age groups first
forvalues s=1/$stageN {
	
	frame change rs_fpm_`s'
	
	forvalues i=1/5 {
		
		preserve
		keep timevar rs_fpm_m_a`i'* rs_fpm_f_a`i'*
		frame copy rs_fpm_`s' rs_fpm_`s'_a`i', replace
		restore
		
	}
	
}

// one stage/age at a time
forvalues s=1/$stageN {
	
	// rubin's rules
	forvalues i=1/5 {
		
		frame change rs_fpm_`s'_a`i'
		
		forvalues k=1/$nimp {
			
			gen ln_rs_fpm_m_a`i'_s`s'_`k' = ln(rs_fpm_m_a`i'_s`s'_`k')
			gen ln_rs_fpm_f_a`i'_s`s'_`k' = ln(rs_fpm_f_a`i'_s`s'_`k')
			
			// Approximate variance on log scale.
			gen var_ln_rs_fpm_m_a`i'_s`s'_`k' = ///
				((ln(rs_fpm_m_a`i'_s`s'_`k'_uci)-ln(rs_fpm_m_a`i'_s`s'_`k'))/1.96)^2
			gen var_ln_rs_fpm_f_a`i'_s`s'_`k' = ///
				((ln(rs_fpm_f_a`i'_s`s'_`k'_uci)-ln(rs_fpm_f_a`i'_s`s'_`k'))/1.96)^2
			
		}
	
		// take the mean and exponentiate
		egen ln_rs_fpm_m_a`i'_s`s' = rowmean(ln_rs_fpm_m_a`i'_s`s'*)
		gen rs_fpm_m_a`i'_s`s' = exp(ln_rs_fpm_m_a`i'_s`s')
		
		egen ln_rs_fpm_f_a`i'_s`s' = rowmean(ln_rs_fpm_f_a`i'_s`s'*)
		gen rs_fpm_f_a`i'_s`s' = exp(ln_rs_fpm_f_a`i'_s`s')
		
		// calculate se and CIs
		egen ln_rs_fpm_m_a`i'_s`s'_B = rowsd(ln_rs_fpm_m_a`i'_s`s'_*) 
		egen ln_rs_fpm_m_a`i'_s`s'_U = rowmean(var_ln_rs_fpm_m_a`i'_s`s'_*)
		gen ln_rs_fpm_m_a`i'_s`s'_se = ///
			sqrt(ln_rs_fpm_m_a`i'_s`s'_U+(1+1/${nimp})*ln_rs_fpm_m_a`i'_s`s'_B^2)
		gen rs_fpm_m_a`i'_s`s'_lci = ///
			exp(ln_rs_fpm_m_a`i'_s`s'-1.96*ln_rs_fpm_m_a`i'_s`s'_se)
		gen rs_fpm_m_a`i'_s`s'_uci = ///
			exp(ln_rs_fpm_m_a`i'_s`s'+1.96*ln_rs_fpm_m_a`i'_s`s'_se)
		
		egen ln_rs_fpm_f_a`i'_s`s'_B = rowsd(ln_rs_fpm_f_a`i'_s`s'_*) 
		egen ln_rs_fpm_f_a`i'_s`s'_U = rowmean(var_ln_rs_fpm_f_a`i'_s`s'_*)
		gen ln_rs_fpm_f_a`i'_s`s'_se = ///
			sqrt(ln_rs_fpm_f_a`i'_s`s'_U+(1+1/${nimp})*ln_rs_fpm_f_a`i'_s`s'_B^2)
		gen rs_fpm_f_a`i'_s`s'_lci = ///
			exp(ln_rs_fpm_f_a`i'_s`s'-1.96*ln_rs_fpm_f_a`i'_s`s'_se)
		gen rs_fpm_f_a`i'_s`s'_uci = ///
			exp(ln_rs_fpm_f_a`i'_s`s'+1.96*ln_rs_fpm_f_a`i'_s`s'_se)
		
	}
	
}

// delete any variables not needed
// females
forvalues s=1/$stageN {
	
	// rubin's rules
	forvalues i=1/5 {
		
		frame change rs_fpm_`s'_a`i'
		
		drop ln_rs_fpm_f_a`i'_s`s' ln_rs_fpm_f_a`i'_s`s'_B ln_rs_fpm_f_a`i'_s`s'_U ln_rs_fpm_f_a`i'_s`s'_se
		
		forvalues k=1/$nimp {
		
			drop rs_fpm_f_a`i'_s`s'_`k' rs_fpm_f_a`i'_s`s'_`k'_uci rs_fpm_f_a`i'_s`s'_`k'_lci ln_rs_fpm_f_a`i'_s`s'_`k' var_ln_rs_fpm_f_a`i'_s`s'_`k'
			
		}
			
	}
	
}

// males
forvalues s=1/$stageN {
	
	forvalues i=1/5 {
		
		frame change rs_fpm_`s'_a`i'
		
		drop ln_rs_fpm_m_a`i'_s`s' ln_rs_fpm_m_a`i'_s`s'_B ln_rs_fpm_m_a`i'_s`s'_U ln_rs_fpm_m_a`i'_s`s'_se
		
		forvalues k=1/$nimp {
		
			drop rs_fpm_m_a`i'_s`s'_`k' rs_fpm_m_a`i'_s`s'_`k'_uci rs_fpm_m_a`i'_s`s'_`k'_lci ln_rs_fpm_m_a`i'_s`s'_`k' var_ln_rs_fpm_m_a`i'_s`s'_`k'
			
		}
			
	}
	
}

// merge the stage and age frames into one
// merge over stage first
forvalues i=1/5 {
	frame copy rs_fpm_1_a`i' rs_fpm_a`i'
}	

forvalues i=1/5 {
	
	frame change rs_fpm_a`i'
		
	forvalues s=2/$stageN {
	
		frlink m:1 timevar, frame(rs_fpm_`s'_a`i')
			
		frget *, from(rs_fpm_`s'_a`i') exclude(timevar)
	
	}
	
}

// then merge over age
frame copy rs_fpm_a1 rs_fpm
frame change rs_fpm

forvalues i=2/5 {
	
	frlink m:1 timevar, frame(rs_fpm_a`i')
			
	frget *, from(rs_fpm_a`i') exclude(timevar)
		
}

save "rs_fpm", replace

frame change default

*** PLOTS ***
*** CIs - one row for each sex ***

local ytitle_m_1 "Male"
local ytitle_f_1 "Female"
local yaxistitle_1 "Relative Survival"

forvalues i=1/5 {
	
	// male
	
	frame rs_fpm:                                              	 			 ///
		twoway (rarea rs_fpm_m_a`i'_s1_lci rs_fpm_m_a`i'_s1_uci 			 ///
					  time													 ///
					  	  , 										   		 ///
					  	  sort												 ///
					  	  color(lavender%50)	 							 ///
					  	  lwidth(none) 										 ///
					  	  fintensity(inten20)) 								 ///
				(rarea rs_fpm_m_a`i'_s2_lci rs_fpm_m_a`i'_s2_uci 			 ///
					  time													 ///
					  	  , 										   		 ///
					  	  sort												 ///
					  	  color("173 40 113"%50)	 						 ///
					  	  lwidth(none) 										 ///
					  	  fintensity(inten20))  							 ///
				(rarea rs_fpm_m_a`i'_s3_lci rs_fpm_m_a`i'_s3_uci 			 ///
					  time													 ///
					  	  , 										   		 ///
					  	  sort												 ///
					  	  color(navy%50)	 							 	 ///
					  	  lwidth(none) 										 ///
					  	  fintensity(inten20)) 								 ///
			(line rs_fpm_m_a`i'_s1 rs_fpm_m_a`i'_s2                  		 ///
					rs_fpm_m_a`i'_s3 				                  		 ///
					time                                               		 ///
						,                                              		 ///
						sort                                           		 ///
						lcolor(lavender "173 40 113" navy))   	 			 ///
							, 												 ///
							name(age`i'_ci_m, replace)                       ///
							aspectratio(1)                              	 ///
							xtitle(Time after diagnosis (years))        	 ///
							xscale(range(0 10))                         	 ///
							xlabel(0(2)10)                              	 ///
							l1title("`ytitle_m_`i'' ", size(medlarge)) 		 ///
							ytitle("`yaxistitle_`i'' ")                   	 ///
							yscale(range(0 1.6))                          	 ///
							ylabel(0(0.2)1, format(%2.1f))              	 ///
							scale(0.8)
	
	frame rs_pp: addplot:  											 		 ///
			(rcap rs_m_a`i'_s1_lci rs_m_a`i'_s1_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color(lavender))					    				 ///
			(rcap rs_m_a`i'_s2_lci rs_m_a`i'_s2_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color("173 40 113"))					    			 ///
			(rcap rs_m_a`i'_s3_lci rs_m_a`i'_s3_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color(navy))					    			 		 ///
			(scatter rs_m_a`i'_s1 rs_m_a`i'_s2                               ///
					 rs_m_a`i'_s3                          					 ///
					 time                                                    ///
						,                                                    ///
						pstyle(p1line p2line p3line)                 		 ///
						mcolor(lavender "173 40 113" navy))     			 ///
							,												 ///
							yscale(range(0 1.6))                             ///
							ylabel(0(0.2)1, format(%2.1f))                   ///
							legend(order(4 "FPM - SI" 5 "FPM - SII" 		 ///
									 6 "FPM - SIII"   			 			 ///
									 10 "PP - SI" 11 "PP - SII"    			 ///
									 12 "PP - SIII") 	 		 		     ///
								row(2) pos(6) ring(3) size(vsmall)         	 ///
									symxsize(*0.5))
			
	// female
			
	frame rs_fpm:                                            		 		 ///
		twoway (rarea rs_fpm_f_a`i'_s1_lci rs_fpm_f_a`i'_s1_uci	 	 		 ///
			  time													 		 ///
			  	  , 										   		 		 ///
			  	  sort												 		 ///
			  	  color(lavender%50)	 						 	 		 ///
			  	  lwidth(none) 										 		 ///
			  	  fintensity(inten20)) 								 		 ///
		(rarea rs_fpm_f_a`i'_s2_lci rs_fpm_f_a`i'_s2_uci 			 		 ///
			  time													 		 ///
			  	  , 										   		 		 ///
			  	  sort												 		 ///
			  	  color("173 40 113"%50)						     		 ///
			  	  lwidth(none) 										 		 ///
			  	  fintensity(inten20)) 								 		 ///
		(rarea rs_fpm_f_a`i'_s3_lci rs_fpm_f_a`i'_s3_uci 			 		 ///
			  time													 		 ///
			  	  , 										   		 		 ///
			  	  sort												 		 ///
			  	  color(navy%50)	 								 		 ///
			  	  lwidth(none) 										 		 ///
			  	  fintensity(inten20)) 								 		 ///  
		(line rs_fpm_f_a`i'_s1 rs_fpm_f_a`i'_s2                  	 		 ///
			  rs_fpm_f_a`i'_s3                	 		 					 ///
			  time                                               	 		 ///
					,                                              	 		 ///
					sort                                           	 		 ///
					lcolor(lavender "173 40 113" navy)) 		 			 ///
						, 											 		 ///
						name(age`i'_ci_f, replace)                   		 ///
						aspectratio(1)                               		 ///
						xtitle(Time after diagnosis (years))         		 ///
						xscale(range(0 10))                          		 ///
						xlabel(0(2)10)                               		 ///
						l1title("`ytitle_f_`i'' ", size(medlarge)) 	 		 ///
						ytitle("`yaxistitle_`i'' ")                  		 ///
						yscale(range(0 1.6))                           		 ///
						ylabel(0(0.2)1, format(%2.1f))               		 ///
						scale(0.8)
						
		frame rs_pp: addplot:  												 ///
			(rcap rs_f_a`i'_s1_lci rs_f_a`i'_s1_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color(lavender))					    			 	 ///
			(rcap rs_f_a`i'_s2_lci rs_f_a`i'_s2_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color("173 40 113"))					    			 ///
			(rcap rs_f_a`i'_s3_lci rs_f_a`i'_s3_uci time					 ///
					,														 ///
					msize(vsmall)											 ///
					lwidth(thin)											 ///
					color(navy))					    				 	 ///
			(scatter rs_f_a`i'_s1 rs_f_a`i'_s2                               ///
					 rs_f_a`i'_s3                                			 ///
					 time                                                    ///
						,                                                    ///
						pstyle(p1line p2line p3line)                  		 ///
						mcolor(lavender "173 40 113" navy))    		 		 ///
							,												 ///
							yscale(range(0 1.6))                             ///
							ylabel(0(0.2)1, format(%2.1f))                   ///
							legend(order(4 "FPM - SI" 5 "FPM - SII" 		 ///
									 6 "FPM - SIII"   			 			 ///
									 10 "PP - SI" 11 "PP - SII"    			 ///
									 12 "PP - SIII") 	 		 			 ///
								row(2) pos(6) ring(3) size(vsmall)         	 ///
									symxsize(*0.5))
	
}

// ideally - increase y axis for all plots to make scale the same where rs>1

local graph_label_1 Age 18-44 years
local graph_label_2 Age 45-54 years
local graph_label_3 Age 55-64 years
local graph_label_4 Age 65-74 years
local graph_label_5 Age 75+ years

forvalues i=1/5 {
	grc1leg age`i'_ci_m age`i'_ci_f, ///
		col(1) pos(6) iscale(0.8) name(age`i'_ci_2, replace) ///
		imargin(0 0 0 0) title("          `graph_label_`i''", size(medsmall)) ycommon
}

grc1leg age1_ci_2 age2_ci_2 age3_ci_2 age4_ci_2 age5_ci_2, ///
	row(1) pos(6) iscale(0.8) imargin(0 0 0 0) ycommon

graph export "valid_ci.pdf", replace
graph export "valid_ci.png", replace













