* Reproduces Tables 3 and 4.
* Input: data_incidence.csv.
* Outputs: Table3.tex and Table4.tex.

insheet using data_incidence.csv, clear name
rename city_d city

* Subsidy Distributional Effect
gen ln_survey_income = ln(survey_income + 0.001)
gen ev_income = ev*survey_income
gen ev_lnincome = ev*ln_survey_income

gen survey_inc_k = survey_income/10^3
gen s_c_k = s_c*10

* no FE
reghdfe s_c_k ln_survey_income if ev, vce(r)
est store m1
estadd local CityFE "No"
estadd local TimeFE "No"
estadd local InterFE "No"

* city FE
reghdfe s_c_k ln_survey_income if ev, absorb(city) vce(r)
est store m2
estadd local CityFE "Yes"
estadd local TimeFE "No"
estadd local InterFE "No"

* hy FE
reghdfe s_c_k ln_survey_income if ev, absorb(hy) vce(r)
est store m3
estadd local CityFE "No"
estadd local TimeFE "Yes"
estadd local InterFE "No"

* city * hy FE
reghdfe s_c_k ln_survey_income if ev, absorb(city#hy) vce(r)
est store m4
estadd local CityFE "Yes"
estadd local TimeFE "Yes"
estadd local InterFE "Yes"

esttab m1 m2 m3 m4 ///
    using "Table4.tex", ///
    replace booktabs ///
    se b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(CityFE TimeFE InterFE N r2, ///
          fmt(%9s %9s %9s %9.0fc %9.3f) ///
          labels("City FE" "Time FE" "City×Time FE" "Observations" "R-squared")) ///
    mtitle("(1)" "(2)" "(3)" "(4)") ///
    label ///
    nonumbers ///
    prehead("\begin{table}[htbp]" ///
            "\centering" ///
            "\caption{Stylized Facts about the Correlation between Subsidies and Individual Incomes}" ///
            "\begin{tabular}{lcccc}" ///
            "\toprule" ///
            "& (1) & (2) & (3) & (4) \\" ///
            "\midrule") ///
    postfoot("\bottomrule" ///
             "\end{tabular}" ///
             "\end{table}")
	
* summary statistics ******************************
insheet using data_incidence.csv, clear name
eststo clear
gen subsidy_all_k = (s_c + subsidy_local_net)*10
gen msrp_k = msrp*10
* summary
foreach var in sales msrp_k power cost_unit weight size subsidy_all_k import at ev suv {
    qui sum `var'
    matrix stats = nullmat(stats) \ ///
        [`r(N)', `r(mean)', `r(sd)', `r(min)', `r(max)']
    local varlabels `"`varlabels' "`:var label `var''""'
}

//esttab matrix(stats, fmt(%12.0fc %12.2f %12.2f %12.2f %12.2f)) 
* generate table
esttab matrix(stats, fmt(%12.0fc %12.2f %12.2f %12.2f %12.2f)) ///
    using "Table3.tex", ///
    mlabels("Observations" "Mean" "Std. Dev." "Min" "Max") ///
    coeflabels(`varlabels') ///
    booktabs replace


