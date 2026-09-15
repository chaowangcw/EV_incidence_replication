* Reproduces the second column of Table A1 without PHEVs.
* Input: data_incidence.csv.
* Output: TableA1col2.rtf.

insheet using data_incidence.csv, clear name

preserve

* Calculate city-level averages first
collapse (mean) avg_price=msrp (first) p_id [aw=sales], by(city_d ev)

* Create vehicle type variable
gen vehicle_type = "ICEV" if ev == 0
replace vehicle_type = "EV" if ev == 1
//replace vehicle_type = "PHEV" if phev == 1

* Calculate summary statistics by vehicle type
eststo clear
estpost tabstat avg_price, by(vehicle_type) statistics(mean sd) columns(statistics)

* Store the city and province counts from original data
restore
distinct city_d
local city_count = r(ndistinct)
distinct p_id
local province_count = r(ndistinct)

* Calculate shares in original data
tempvar total_sales ev_sales
egen `total_sales' = total(sales)
egen `ev_sales' = total(sales) if ev == 1
//egen `phev_sales' = total(sales) if phev == 1

summ `total_sales'
local total = r(mean)
summ `ev_sales'
local ev_share = r(mean)/`total'
//summ `phev_sales'
//local phev_share = r(mean)/`total'
local icev_share = 1 - `ev_share'

* Create a matrix for the additional statistics
matrix additional = J(8, 2, .)
matrix rownames additional = "Cities" "Provinces" "EV_Share" "ICEV_Share" "EV_Price" "ICEV_Price" "EV_SD" "ICEV_SD"
matrix colnames additional = "Statistic" "Value"

matrix additional[1,2] = `city_count'
matrix additional[2,2] = `province_count'
matrix additional[3,2] = `ev_share'
//matrix additional[4,2] = `phev_share'
matrix additional[4,2] = `icev_share'

* Get the price statistics from estpost results
matrix price_stats = e(mean)
matrix list price_stats
matrix sd_stats = e(sd)
matrix list sd_stats

matrix additional[5,2] = price_stats[1,1]  // EV mean
//matrix additional[7,2] = price_stats[1,2]  // PHEV mean  
matrix additional[6,2] = price_stats[1,2]  // ICEV mean
matrix additional[7,2] = sd_stats[1,1]     // EV SD
//matrix additional[10,2] = sd_stats[1,2]    // PHEV SD
matrix additional[8,2] = sd_stats[1,2]    // ICEV SD

* Export to Word using esttab
esttab matrix(additional) using "TableA1col2.rtf", replace ///
    title("Summary Statistics") ///
    cells(b(fmt(%9.0g %9.4f %9.2f))) ///
    noobs note("Notes: Price statistics are averages of city-level average prices.")
