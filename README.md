# CMSY Multi Species forecast
Multi species forecasting method based on the CMSY stock assessment method's output. The CSV files report examples of inputs to the method, e.g. for multiple species (Out_March062019_O_Stocks_ID_19_Med.csv), single species (Out_March062019_O_Stocks_ID_19_Med_single.csv and Out_SingleSpecies.csv), aggregated groups of species (Out_Ukr_ID.csv).

## CMSY_recovery_forecasts_legacy.R ##
Source code used for the experiments reported in the paper
["Froese, R., Winker, H., Coro, G., Demirel, N., Tsikliras, A. C., Dimarchopoulou, D., ... & Matz-Lück, N. (2018). Status and rebuilding of European fisheries. Marine Policy, 93, 159-170."](https://www.sciencedirect.com/science/article/pii/S0308597X17307364?via%3Dihub)

which assumes that F/Fmsy can be set for all species to a fixed value (e.g. 0.5, 0.8, 0.95) for the forecast years.

## CMSY_recovery_forecasts_based_on_last_year_f_fmsy.R ##
Code used for the experiments reported in the paper ["Armelloni, E. N., Scanu, M., Masnadi, F., Coro, G., Angelini, S., and Scarcella, G. 2021. Data Poor Approach for the Assessment of the Main Target Species of Rapido Trawl Fishery in Adriatic Sea. Frontiers in Marine Science, 8."](https://www.frontiersin.org/articles/10.3389/fmars.2021.552076/full)

which assumes that F/Fmsy can be set for all species to a portion of the value in the latest year (e.g. 0.5 F_latest_year/Fmsy, etc.) for the forecast years.

## CMSY_recovery_short_term_forecasts.R ##

This version builds on "CMSY_recovery_forecasts_based_on_last_year_f_fmsy.R". It accept short term forecasts and single species dataset (i.e. Out_SingleSpecies.csv). In the case of short term forecast (from 1 to 5 years), scenarios are not activated and projections base on Fcur for last year.
