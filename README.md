**df_HAB.R:** sample wind speed data from WRF output for Hafar-Al-Batin, Saudi Arabia, for testing functionality of extrapolation function

**windextrap_harmonics.R:** function fitting wind extrapolation model, manual below.

Created by M.Alifa

Final version for publication, January 7 2021

This function fits a nonlinear, temporally-varying model for vertical extrapolation of hourly wind speed for the purpose of wind energy applications, and performs wind speed extrapolation in the presence of the necessary data

More details about the model can be found in Crippa et al, 2021 (https://doi.org/10.1016/j.apenergy.2021.117378)
  	
Inputs:

* "train" is the dataframe with data for model training, with column names:
  * "hour" for the hour of day for each observation
  * "v1"  for time series of near ground wind speeds
  * "h1" for height above ground of each v1 observation
  * "v2"  for time series of wind speeds aloft
  * "h2" for height above ground of each v2 observation
* "test" (optional), a dataframe with the same column names as "train" ("v2" is not necessary)
if specified, the function will return extrapolated wind speeds and their 95% Confidence interval for the testing set data

Output:

A list with the following variables:
* "model": class 'nls', harmonics model for wind speed extrapolation
* "res_model": class 'nls', harmonics model for the error variance

Additional output when "test" is specified:
* "v2pred": vector of extrapolated wind speeds
* "v2_lower95": vector of lower 95% CI for the extrapolated wind speed
* "v2_upper95": vector of upper 95% CI for the extrapolated wind speed

Function requires the following R libraries: nlstools

