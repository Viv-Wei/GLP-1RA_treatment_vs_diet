# Raw data

The files below are not redistributed. Download them and place them at the paths shown before running the `params` or `main` stages. 

## Global Burden of Disease Study 2023 (IHME)

GBD Results Tool, https://vizhub.healthdata.org/gbd-results/ (free registration; IHME terms of use apply).
Common settings: year 2023; sex male, female, both; age 25–29 to 95+ years in 5-year groups; locations global and
Australia, Bangladesh, Brazil, Bulgaria, Canada, China, Egypt, Germany, India, Israel, Japan, Mexico, New Zealand,
Nigeria, Romania, South Africa, Türkiye, United Kingdom, United States of America.

| File | Measure | Cause | Risk factor | Metric |
|---|---|---|---|---|
| `data/gbd/queryA_ckd_3years.csv` | Deaths, DALYs | Chronic kidney disease | Dietary risks and its 13 component dietary risks | Number, percent, rate |
| `data/gbd/queryA_ICH_all_merged.csv` | Deaths, DALYs | Intracerebral hemorrhage | As above | Number, percent, rate |
| `data/gbd/queryA_SAH_all_merged.csv` | Deaths, DALYs | Subarachnoid hemorrhage | As above | Number, percent, rate |
| `data/gbd/cause_totals/queryB_t2d.csv` | Deaths, incidence, prevalence | Diabetes mellitus type 2 | None (all causes of the outcome) | Number, rate |
| `data/gbd/1-3_population.csv` | Population | None | None | Number |

GBD 2023 risk factor exposure estimates, Global Health Data Exchange (https://ghdx.healthdata.org/):

| File | Dataset |
|---|---|
| `data/gbd/bmi_exposure/IHME_GBD_2023_RISK_EXPOSURE_HIGH_BMI_PREVALENCE_OF_OBESITY_AND_OVERWEIGHT_Y2025M10D10.CSV` | GBD 2023 risk exposure, high BMI, prevalence of overweight and obesity |

## UN World Population Prospects 2024

`data/wpp/WPP2024_PopulationByAge5GroupSex_Medium.csv`, population by 5-year age group and sex, medium variant,
https://population.un.org/wpp/downloads (CSV format). Used for the population-weighted global overweight prevalence
and the 2053 population sensitivity analysis.

## Wittgenstein Centre Human Capital Data Explorer

`data/wittgenstein/wcde_clean.csv`, population by age, sex, and educational attainment, 2020, all countries,
https://dataexplorer.wittgensteincentre.org/ (columns Area, Year, Age, Sex, Education, Population). Used for the
education structure of adults aged 25 years and older.

## Global Dietary Database

`data/exposure_19country.csv`, dietary intake by country, age, sex, and education level for the 19 study countries,
https://www.globaldietarydatabase.org/ (data access subject to GDD terms). Used for education-specific
diet-attributable fractions.
