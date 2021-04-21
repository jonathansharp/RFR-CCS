% This script runs a number of functions to create a multi-year monthly
% surface pCO2 climatology in a region that includes the California
% Current. Collectively, the scrpit (1) imports SOCATv2020 observations of
% surface pCO2, (2) grids those observations to a monthly 0.25x0.25 degree
% grid, (3) pulls in satellite and model data (SST, SSS, CHL, wind speed,
% etc.) to co-locate with pCO2 grid cells, and (4) uses those co-located
% variables to interpolate surface surface pCO2 to a continuous 0.25x0.25
% degree map from 1998 to 2019.

% this can be changed to split training/test data randomly (1), by
% withholding specific moorings (2), or by withholding specific years (3)
split = 0;
startyear = [1995,1996,1997,1998,1999];
% this can be changed to obtain different random splits of training/test data
rng_seed = 8;
numsplits = 20;
% moorings to withhold
omitmoors = {'WA' 'CCE1' 'CCE2' 'SEAK' 'NH10' 'CB06' 'LaPush' 'KwakshuaChannel' 'Dabob' 'Exp'};
% this can be changed to make predictions on grid (1) or not (0)
grid = 1;
grid_val = 1;
% som_map = 0;

% This imports and splits data into training/test
importSOCATv2020;
% This grids SOCAT data
gridSOCATv2020;
% This detrends observations
detrendSOCATv2020;
% This loads predictor variables
loadvarsSOCATv2020;
% This interpolates pCO2 by random forest using training data and tests the
% results against test data
%interpolateSOCATv2020_RF_validate;
% This interpolates pCO2 by random forest using all data
interpolateSOCATv2020_RF_all;
% This interpolates pCO2 by SOM-FFN using training data and tests the
% results against test data
% interpolateSOCATv2020_FFN_validate;
% This interpolates pCO2 by SOM-FFN using all data
% interpolateSOCATv2020_FFN_all;

importLAND
importLAR

% This determines CO2 flux
calculateCO2flux
%calculateCO2flux_moorings


% Percentage of filled grid cells
% Index = ~isnan(SOCATv2020_grid.SSS) & ...
%         ~isnan(SOCATv2020_grid.MLD) & ...
%         ~isnan(SOCATv2020_grid.SST);
% sum(~isnan(SOCATv2020_grid.all.pco2_ave_weighted(Index)))./...
% sum(~isnan(SOCATv2020_grid.SST(Index)))