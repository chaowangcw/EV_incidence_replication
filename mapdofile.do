* Reproduces Figure 2: sales-weighted average local EV subsidies.
* Runtime inputs: china_label.dta (map attributes and labels) and
*                 china_map.dta (polygon coordinates).
* Output: localsubsidy.eps in the replication root.
* Run this do-file with the replication root as Stata's working directory.

clear all
set more off

* Load map attributes and label locations.
use "china_label.dta", clear

* Remove unmatched records and verify one observation per polygon ID.
drop if missing(id)
isid id

spmap x using "china_map.dta", id(id) ///
    label( ///
        label(ename) ///
        xcoord(x_coord) ///
        ycoord(y_coord) ///
        size(*.5) ///
    ) ///
    clmethod(custom) ///
    clbreaks(0 3 6 9 12 99999) ///
    fcolor(Greens2) ///
    ocolor(white ..) ///
    osize(medthin ..) ///
    legend( ///
        size(*1.3) ///
        order(2 "<=3" 3 "3~6" 4 "6~9" 5 "9~12" 6 ">=12") ///
        title("Sales-Weighted Average" "Subsidy (RMB '000)", size(*0.5)) ///
    )

graph export "localsubsidy.eps", as(eps) preview(off) replace
