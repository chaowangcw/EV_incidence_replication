* Reproduces Table A2.
* Input: data_incidence.csv.
* Output: TableA2.tex.

insheet using data_incidence.csv, clear name
rename city_d city

bysort model_variant hy city: egen avg_subsidy = mean(s_c)
gen ln_subsidy = ln(avg_subsidy+ .001)
gen volume = length*width*height/10^9

preserve 

* Create unique observations by model_variant, hy, and city
bysort model_variant hy city: gen tem1 = _n == 1
keep if tem1 

* Run regressions and store results
reghdfe msrp avg_subsidy power kw_100km weight volume suv if ev, absorb(brand_d hy) vce(r)
est store model1

reghdfe msrp avg_subsidy power kw_100km weight volume suv hhinc_mu if ev, absorb(brand_d hy) vce(r)
est store model2

reghdfe msrp avg_subsidy power kw_100km weight volume suv hhinc_mu subsidy_local_net if ev, absorb(brand_d hy) vce(r)
est store model3

* Display results side by side
esttab model1 model2 model3, ///
       t stats(r2 N, labels("R-squared" "Observations")) ///
       star(* 0.10 ** 0.05 *** 0.01) ///
       mtitle("Model 1" "Model 2" "Model 3") ///
       label

* Export to LaTeX
esttab model1 model2 model3 using "TableA2.tex", ///
       se ar2 stats(r2 N, labels("R-squared" "Observations")) ///
       star(* 0.10 ** 0.05 *** 0.01) ///
       mtitle("Model 1" "Model 2" "Model 3") ///
       label replace

restore

