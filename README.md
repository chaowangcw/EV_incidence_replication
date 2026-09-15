# Replication Package: The Incidence and Distributional Effects of Electric Vehicle Subsidies in China

## 1. Overview

This directory contains the MATLAB and Stata programs together with an approximately 20% public-use subsample of the principal analysis data. The subsample covers product-city-half-year observations from 2016 through 2019 across all 13 cities and all 104 city-half-year markets. It is intended to permit inspection and testing of the replication code; it is not the full analysis sample used for the published estimates.

Results obtained from this directory will generally differ from the tables and figures reported in the paper. The published results were estimated using the full 89,104-observation analysis sample, whereas the included `data_incidence.csv` contains 17,831 observations (20.011% of the full row count). The public-use subsample is therefore not intended to provide numerical reproduction of the published estimates.

The main structural specification is implemented in `main_UniformMSRP.m`. Two alternative demand specifications are implemented in `demand_BLPIV.m` and `demand_priceoverincome/demand_priceoverincome.m`. The Stata programs produce the descriptive and reduced-form tables. `mapdofile.do` produces Figure 2, and `graph_rep.m` contains the retained plotting code for Figures 1, 3, and 4.

All text files should be opened as UTF-8.

## 2. Citation and package version

### Paper citation

> Jing Liang, Chao Wang, and Junji Xiao. "The Incidence and Distributional Effects of Electric Vehicle Subsidies in China" *[Journal of Economic Behavior & Organization]*.

### Package version

| Field | Entry |
|---|---|
| Release date | September 2026 |

## 3. Directory structure

```text
replication/
|-- README.md
|-- data_incidence.csv
|-- run_all_matlab.m
|-- main_UniformMSRP.m
|-- demand_BLPIV.m
|-- graph_rep.m
|-- mapdofile.do
|-- china_label.dta
|-- china_map.dta
|-- NEVproduction.xls
|-- subsidy&income.xls
|-- Table3and4.do
|-- TableA1.do
|-- TableA1col2.do
|-- TableA2.do
|-- econometric toolbox/
|   |-- distrib/beta_inv.m
|   |-- distrib/beta_pdf.m
|   |-- distrib/tdis_inv.m
|   |-- regress/ols.m
|   |-- regress/tsls.m
|   `-- util/cols.m
|-- matlab2tikz-matlab2tikz-806c97d/
|   |-- LICENSE.md
|   `-- src/
|-- demand_priceoverincome/
|   |-- demand_priceoverincome.m
|   `-- specification-specific MATLAB functions
|-- diary/
|   `-- MATLAB text and MAT outputs
`-- other MATLAB function files
```

The MATLAB function files in the root directory support the main specification and the BLP-IV alternative. The same-named functions under `demand_priceoverincome/` implement the price-over-income specification and intentionally take precedence when that script is run.

## 4. Data files

### Included data

`data_incidence.csv` is the approximately 20% public-use subsample used by the main MATLAB programs and by `Table3and4.do`, `TableA1col2.do`, and `TableA2.do` in this test directory. It contains 17,831 observations for 2016-2019. The unit of observation is a vehicle product in a city and half-year market.

### Construction and limitations of the 20% subsample

The release subsample was constructed by selecting `model_variant` values globally and retaining all available city and half-year observations for each selected model. Thus, a selected model is not independently included or excluded city by city. Models needed by the Figure 5 elasticity calculation were included with certainty, and the remaining models were selected using random seed `20260903` until the retained row count reached approximately 20% of the full data.

Because the Figure 5 models are concentrated in later periods, the sampling fraction is not exactly 20% within each half-year. The retained row shares range from 16.71% in `hy=1` to 28.69% in `hy=8`. This feature preserves the required Figure 5 products but makes the public-use file a testing sample rather than a representative 20% sample within every period.

Several variables in `data_incidence.csv` were constructed before subsampling. In particular, the stored instrument columns `mktn` and `iv_local_msrp_hat` (the local-price IV) retain their values from the full analysis data for the selected rows; they are not recomputed using the reduced product choice sets. The same qualification applies to other preconstructed BLP and market/local IV columns stored in the CSV. By contrast, the differentiation instruments created by `build_differentiation_iv.m` are recalculated by MATLAB from the included subsample.

The stored product market shares and other row-level constructed variables likewise retain their full-data values. Product shares are not renormalized after unselected products are removed. When the demand programs calculate the outside option from the retained products, omitted products consequently become part of the residual outside share. This changes the effective choice set relative to the full-sample estimation.

Subsampling also creates sparsely supported fixed effects. The demand drivers therefore omit brand dummy columns supported by fewer than two observations and remove micro-instrument columns with fewer than two nonzero observations. These are numerical-stability adjustments for the public-use subsample and are another reason its estimates need not equal the published full-sample estimates.

`NEVproduction.xls` supplies the monthly EV and total vehicle production series used for Figure 1. `subsidy&income.xls` supplies the income-group and subsidy data used for Figure 3; the retained specification in `graph_rep.m` reads `sheet1`.

Figure 2 uses two included Stata files. `china_label.dta` contains the sales-weighted average local EV subsidy variable `x`, polygon identifiers, English labels, and label coordinates. `china_map.dta` contains the corresponding polygon coordinates.

### Data availability and provenance statement

This test package contains a public-use subsample of the principal analysis file and the ancillary files needed to exercise the code, except for the restricted full-sample file used by `TableA1.do`. It does not contain the full data needed for numerical reproduction of the reported results. The authors constructed the analysis files from the vehicle-sales, product-characteristics, income, subsidy-policy, production, and geographic sources described below. The raw private records are not included.

| Replication file | Provenance and construction | Availability and use |
|---|---|---|
| `data_incidence.csv` | This author-constructed public-use file is an approximately 20% model-based subsample of the 89,104-observation analysis data. It contains 17,831 observations, 13 cities, and eight half-year periods from 2016 through 2019. The underlying full analysis file combines vehicle sales, product characteristics, subsidy schedules, city-market variables, and income measures. Vehicle sales are measured using compulsory vehicle-insurance policy records from the private data source identified in the paper. The income measures combine resident-income statistics with the passenger-vehicle buyer survey conducted by CVSC-TNS Research (CTR), as described in the paper. Constructed variables stored in the released rows generally retain their full-data values. | Included for code inspection and testing. It is the common test input for the three MATLAB demand specifications and for `Table3and4.do`, `TableA1col2.do`, and `TableA2.do`. Estimates obtained from it are not expected to equal the published full-sample results. |
| `NEVproduction.xls` | This workbook contains the monthly EV and total vehicle production series used to construct Figure 1. | Included and read directly by the Figure 1 section of `graph_rep.m`. |
| `subsidy&income.xls` | This author-constructed workbook contains the income-group and subsidy values used for Figure 3. The retained replication code reads `sheet1`. The values summarize the relationship between consumer income and EV subsidies using the income and subsidy sources described in the paper. | Included and read directly by the Figure 3 section of `graph_rep.m`. |
| `china_label.dta` and `china_map.dta` | These files form the map-ready input for Figure 2. `china_map.dta` stores polygon coordinates. `china_label.dta` stores matching polygon identifiers, English labels, label coordinates, and the author-constructed sales-weighted average local EV subsidy variable `x`. The subsidy values are derived from the local subsidy information used in the analysis. | Included and read by `mapdofile.do`. |

### Data not included

`TableA1.do` requires the restricted file `data_hy_city_no_miss_pct_full_2025.dta`, which is not included in this replication package. It is constructed from the private vehicle-sales data described in the paper and reproduces the full-sample statistics used in Table A1. The authors are not authorized to redistribute these private record-level data. Researchers who independently obtain access from the data source identified in the paper can supply the full path to the file when running the program; no author-specific path is embedded in the code.

| Restricted-data field | Information |
|---|---|
| Exact filename | `data_hy_city_no_miss_pct_full_2025.dta` |
| Reason redistribution is not permitted | The file is constructed from private record-level vehicle-sales data supplied for research use. The authors do not have permission to redistribute the underlying data. |
| Original data provider | The private vehicle-sales data source identified in the Data section of the paper. Sales are measured using compulsory vehicle-insurance policy records, as described in the paper. |
| How to apply for access | Contact the data source identified in the paper directly and request access to the corresponding vehicle-sales/compulsory-insurance records. |

The authors cannot grant access through this package. Access must be requested from the original provider under its applicable terms. The resulting `TableA1.rtf` contains aggregate statistics; it should be included with the package only if its disclosure is permitted by the data-use agreement.

Figure 4 requires variables, `s_inc` and `demogr`, contained in the dated `diary/results_*.mat` file generated by `main_UniformMSRP.m`.

## 5. Software and computational requirements

### MATLAB

The code requires MATLAB and uses functions from:

- Statistics and Machine Learning Toolbox, including `dataset`, `tabulate`, `norminv`, `prctile`, `grpstats`, and `boxplot`;
- Optimization Toolbox, including `fminunc` and `fmincon`; and
- Parallel Computing Toolbox, including `parpool` and `parfor`.

The package includes the required `ols`, `tsls`, `cols`, `tdis_inv`, `beta_inv`, and `beta_pdf` routines under `econometric toolbox/`; this directory is added automatically by each organized MATLAB driver.

The main scripts request a 24-worker parallel pool and use 500 simulation draws per market. Users with fewer cores should change `parpool(24)` to a feasible number. 

`graph_rep.m` uses the bundled `matlab2tikz` snapshot identified as `v1.1.0-99-g806c97d`. The script automatically adds its `src/` directory to the MATLAB path; users therefore do not need to download it separately.

Third-party software information:

- Official distribution page: [matlab2tikz/matlab2tikz on MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz).
- Upstream source repository: [matlab2tikz/matlab2tikz on GitHub](https://github.com/matlab2tikz/matlab2tikz).
- Suggested citation: Nico Schlömer, *matlab2tikz/matlab2tikz*, MATLAB Central File Exchange, retrieved September 8, 2026.
- License: BSD 2-Clause. The bundled copy retains the original copyright notice, license terms, and disclaimer in `matlab2tikz-matlab2tikz-806c97d/LICENSE.md`.

The bundled `matlab2tikz` source has not been modified by the authors of this replication package.

### Stata

The Stata programs were prepared for and run using Stata/MP 18 for Windows. They require the following user-written commands:

- `reghdfe` and its dependency `ftools`;
- `estout`, which provides `esttab`, `eststo`, `estpost`, and `estadd`; and
- `distinct`; and
- `spmap`, used to produce Figure 2.

They can be installed once from Stata with:

```stata
ssc install ftools
ssc install reghdfe
ssc install estout
ssc install distinct
ssc install spmap
```

## 6. Instructions for replication

Make a fresh copy of the package before running the programs. The MATLAB estimation routines create and overwrite intermediate MAT files in their working directories.

To run the MATLAB programs in sequence, change to this directory and run:

```matlab
run_all_matlab
```

The runner executes `main_UniformMSRP.m`, `demand_BLPIV.m`, `demand_priceoverincome/demand_priceoverincome.m`. Completion means that the code has run on the public-use test sample; it does not imply that the resulting numerical estimates reproduce the paper.

### Preferred MATLAB specification and counterfactuals

Open MATLAB and run the complete script:

```matlab
run('FULL_PATH_TO_REPLICATION/main_UniformMSRP.m')
```

The script identifies its own directory, changes to the replication root, adds the bundled econometric toolbox, fixes the random-number states, estimates demand, recovers markups, runs the subsidy counterfactuals, and produces the formatted table output and figures. Run the sections in order because the program relies on global variables and intermediate MAT files.

### Alternative demand specification: BLP instruments

Run:

```matlab
run('FULL_PATH_TO_REPLICATION/demand_BLPIV.m')
```

This program follows the preferred demand-estimation logic but uses traditional BLP instruments. It produces only the BLP-IV columns needed to supplement the preferred differentiation-IV results in Table A3.

### Alternative demand specification: price over income

Run:

```matlab
run('FULL_PATH_TO_REPLICATION/demand_priceoverincome/demand_priceoverincome.m')
```

The script changes to its own subdirectory so that its specification-specific functions take precedence, while reading the common data and econometric toolbox from the replication root. It produces Table A4.

### Stata tables

Start Stata, set the working directory to the replication root, and run the required programs:

```stata
cd "FULL_PATH_TO_REPLICATION"
do Table3and4.do
do TableA2.do
do TableA1col2.do
```

The restricted-data portion of Table A1 can be run only by authorized users who have obtained the required input. Pass the full input-file path as an argument while the replication root remains the current working directory:

```stata
do TableA1.do "FULL_PATH_TO/data_hy_city_no_miss_pct_full_2025.dta"
```

If the argument is omitted or the file cannot be found, the program stops with an explanatory error message. `TableA1col2.do` uses the included public replication data and does not require the restricted file.

### Figure 2

Start Stata, set the working directory to the replication root, and run:

```stata
cd "FULL_PATH_TO_REPLICATION"
do mapdofile.do
```

The program checks that `spmap`, `china_label.dta`, and `china_map.dta` are available, maps the subsidy variable `x` using the manuscript's five class intervals, and exports `localsubsidy.eps` to the replication root.

### Figures 1, 3, and 4

Open `graph_rep.m` and run its three figure sections. The script changes to the replication directory and adds the bundled `matlab2tikz/src` directory automatically. Before running Figure 4, first run `main_UniformMSRP.m`; `graph_rep.m` automatically selects the most recent `diary/results_*.mat` file.

- the EV-production section exports `evproduction.tex` for Figure 1;
- the subsidy-by-income section uses sheet1 of the included `subsidy&income.xls` and exports `incomesubsidygraph.tex` for Figure 3; and
- the buyer-versus-population income section exports `incomeandsimulated.tex` for Figure 4.

The Figure 1 and Figure 3 inputs are included in the replication root. Figure 4 uses the main program's saved results rather than a separately supplied legacy MAT file.

## 7. Table and figure crosswalk

| Paper item | Program | Replication output |
|---|---|---|
| Figure 1 | `graph_rep.m` | `evproduction.tex` |
| Figure 2 | `mapdofile.do` | `localsubsidy.eps` |
| Figure 3 | `graph_rep.m` | `incomesubsidygraph.tex` |
| Figure 4 | `graph_rep.m` | `incomeandsimulated.tex` |
| Figure 5 | `main_UniformMSRP.m` | `elasticity.eps` |
| Figure 6 | `main_UniformMSRP.m` | `EVadop_income.eps` |
| Figure 7 | `main_UniformMSRP.m` | `range_ptr_hy.eps` |
| Figure 8 | `main_UniformMSRP.m` | `subsidy_dist.eps` |
| Figure 9 | `main_UniformMSRP.m` | `range_dist.eps` |
| Table 3 | `Table3and4.do` | `Table3.tex` |
| Table 4 | `Table3and4.do` | `Table4.tex` |
| Table 5 | `main_UniformMSRP.m` | formatted MATLAB diary output |
| Table 6 | `main_UniformMSRP.m` | formatted MATLAB diary output |
| Table 7 | `main_UniformMSRP.m` | formatted MATLAB diary output |
| Table 8 | `main_UniformMSRP.m` | formatted MATLAB diary output |
| Table 9 | `main_UniformMSRP.m` | formatted MATLAB diary output |
| Table A1, full sample | `TableA1.do` | `TableA1.rtf` |
| Table A1, subsample | `TableA1col2.do` | `TableA1col2.rtf` |
| Table A2 | `TableA2.do` | `TableA2.tex` |
| Table A3, BLP-IV columns | `demand_BLPIV.m` | timestamped text file under `diary/` |
| Table A4 | `demand_priceoverincome/demand_priceoverincome.m` | timestamped text file under `diary/` |
| Table A5 | `main_UniformMSRP.m` | formatted MATLAB diary output |

## 8. Output locations and reproducibility notes

- Outputs based on `data_incidence.csv` use the 20% public-use subsample and may differ from the paper's full-sample results. Figures 1-3 use separate included inputs; Figures 4-9, and the tables using `data_incidence.csv` may differ.
- Figures 2 and 5-9 are saved as EPS files in the replication root. Figures 1, 3, and 4 are exported as TikZ `.tex` files.
- Tables 3, 4, and A2 are exported as LaTeX files; Table A1 components are exported as RTF files.
- MATLAB estimation logs and MAT results are saved with timestamps in the replication root and/or `diary/`.
- The restricted full-sample portion of Table A1 requires independently obtained data.
