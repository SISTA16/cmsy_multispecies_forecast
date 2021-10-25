# CMSY Multi Species forecast
Multi species forecasting method based on the CMSY stock assessment method's output. The CSV files report examples of inputs to the method, e.g. for multiple species (Out_March062019_O_Stocks_ID_19_Med.csv), single species (Out_March062019_O_Stocks_ID_19_Med_single.csv), aggregated groups of species (Out_Ukr_ID.csv).

## CMSY_recovery_forecasts_legacy.R ##
Source code used for the experiments reported in the paper
"Froese, R., Winker, H., Coro, G., Demirel, N., Tsikliras, A. C., Dimarchopoulou, D., ... & Matz-Lück, N. (2018). Status and rebuilding of European fisheries. Marine Policy, 93, 159-170."

which assumes that F/Fmsy can be set for all species to a fixed value (e.g. 0.5, 0.8, 0.95) for the forecast years.

## CMSY_recovery_forecasts_based_on_last_year_f_fmsy.R ##
Code used for the experiments reported in the paper
"Scanu, M., Armelloni, E. N., Coro, G., Masnadi, F., Angelini, S., & Scarcella, G. (2021). Data poor approach for the assessment of the main target species of rapido trawl fishery in Adriatic Sea. Frontiers in Marine Science, 8, 681."

which assumes that F/Fmsy can be set for all species to a portion of the value in the latest year (e.g. 0.5 F_latest_year/Fmsy, etc.) for the forecast years.
