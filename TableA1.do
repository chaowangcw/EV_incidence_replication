* Reproduces Table A1 using the full sample.
* Restricted input: data_hy_city_no_miss_pct_full_2025.dta.
* The confidential input data are not included and cannot be redistributed by the
* authors. Researchers must request access from the data source identified in
* the paper and comply with the provider's applicable terms.
* Usage: do TableA1.do "FULL_PATH_TO/data_hy_city_no_miss_pct_full_2025.dta"
* Output: TableA1.rtf (written to the current working directory).

// data statistics comparison
args restricted_data_file

if `"`restricted_data_file'"' == "" {
    display as error "TableA1.do requires the full path to the restricted input file."
    display as error `"Usage: do TableA1.do "FULL_PATH_TO/data_hy_city_no_miss_pct_full_2025.dta""'
    exit 601
}

capture confirm file `"`restricted_data_file'"'
if _rc {
    display as error `"Restricted input file not found: `restricted_data_file'"'
    display as error "Required filename: data_hy_city_no_miss_pct_full_2025.dta"
    exit 601
}

use `"`restricted_data_file'"', clear
keep if hy > 8

* Create summary statistics using estpost
preserve

* Calculate city-level averages first
collapse (mean) avg_price=msrp (first) province [aw=sales], by(city2 ev phev)

* Create vehicle type variable
gen vehicle_type = "ICEV" if ev == 0 & phev == 0
replace vehicle_type = "EV" if ev == 1
/* replace vehicle_type = "PHEV" if phev == 1 */

* Calculate summary statistics by vehicle type
eststo clear
estpost tabstat avg_price, by(vehicle_type) statistics(mean sd) columns(statistics)

* Store the city and province counts from original data
restore
distinct city2
local city_count = r(ndistinct)
distinct province
local province_count = r(ndistinct)

* Calculate shares in original data
tempvar total_sales ev_sales phev_sales
egen `total_sales' = total(sales)
egen `ev_sales' = total(sales) if ev == 1
/* egen `phev_sales' = total(sales) if phev == 1 */

summ `total_sales'
local total = r(mean)
summ `ev_sales'
local ev_share = r(mean)/`total'
/* summ `phev_sales' */
/* local phev_share = r(mean)/`total' */
local icev_share = 1 - `ev_share'

* Create a matrix for the additional statistics
matrix additional = J(8, 2, .)
matrix rownames additional = "Cities" "Provinces" "EV_Share" "ICEV_Share" "EV_Price" "ICEV_Price" "EV_SD" "ICEV_SD"
matrix colnames additional = "Statistic" "Value"

matrix additional[1,2] = `city_count'
matrix additional[2,2] = `province_count'
matrix additional[3,2] = `ev_share'
/* matrix additional[4,2] = `phev_share' */
matrix additional[4,2] = `icev_share'

* Get the price statistics from estpost results
matrix price_stats = e(mean)
matrix sd_stats = e(sd)

matrix additional[5,2] = price_stats[1,1]  // EV mean
/* matrix additional[7,2] = price_stats[1,2]  // PHEV mean   */
matrix additional[6,2] = price_stats[1,2]  // ICEV mean
matrix additional[7,2] = sd_stats[1,1]     // EV SD
/* matrix additional[10,2] = sd_stats[1,2]    // PHEV SD */
matrix additional[8,2] = sd_stats[1,2]    // ICEV SD

* Export to Word using esttab
esttab matrix(additional) using "TableA1.rtf", replace ///
    title("Summary Statistics") ///
    cells(b(fmt(%9.0g %9.4f %9.2f))) ///
    noobs note("Notes: Price statistics are averages of city-level average prices.")
