% This script runs a number of functions to create a multi-year monthly
% surface pCO2 climatology in a region that includes the California
% Current. Collectively, the script (1) imports SOCATv2021 observations of
% surface pCO2, (2) grids those observations to a monthly 0.25x0.25 degree
% grid, (3) pulls in satellite and model data (SST, SSS, CHL, wind speed,
% etc.) to co-locate with pCO2 grid cells, and (4) uses those co-located
% variables to interpolate surface surface pCO2 to a continuous 0.25x0.25
% degree map from 1998 to 2019.

%% Optimize model
% rng_seed = randi([1 100],1);
% split = 0;          % split data randomly
% grid_val = 0;       % do not predict pCO2 on grid
% importSOCATv2021;   % import and split data into training/test
% gridSOCATv2021;     % grid training/test data
% detrendSOCATv2021;  % detrend gridded observations
% loadvarsSOCATv2021; % load predictor variables
% interpolateSOCATv2021_RF_optimize;
%                     % interpolate pCO2 by random forest using training data
% clear

%% For Test #1: split training/test data randomly (random 20%)
rng_seed = randi([1 100],1);
split = 1;          % split data randomly
numsplits = 10;     % perform 10 splits
grid_val = 0;       % do not predict pCO2 on grid
importSOCATv2021;   % import and split data into training/test
gridSOCATv2021;     % grid training/test data
detrendSOCATv2021;  % detrend gridded observations
loadvarsSOCATv2021; % load predictor variables
interpolateSOCATv2021_RF_validate;
                    % interpolate pCO2 by random forest using training data

clear

%% For Test #3: split by withholding years (every fifth)
rng_seed = randi([1 100],1);
split = 3;          % split data by withholding years
startyear = [1995,1996,1997,1998,1999];
                    % define start years
grid_val = 0;       % do not predict pCO2 on grid
importSOCATv2021;   % import and split data into training/test
gridSOCATv2021;     % grid training/test data
detrendSOCATv2021;  % detrend gridded observations
loadvarsSOCATv2021; % load predictor variables
interpolateSOCATv2021_RF_validate;
                    % interpolate pCO2 by random forest using training data
clear

%% For Test #2: split by withholding moorings (one by one)
rng_seed = randi([1 100],1);
split = 2;          % split data by withholding moorings
grid_val = 1;       % predict pCO2 on grid
importSOCATv2021;   % import and split data into training/test
gridSOCATv2021;     % grid training/test data
detrendSOCATv2021;  % detrend gridded observations
loadvarsSOCATv2021; % load predictor variables
interpolateSOCATv2021_RF_validate;
                    % interpolate pCO2 by random forest using all data
grid_all = 1;       % predict pCO2 on grid
interpolateSOCATv2021_RF_all;
                    % interpolate pCO2 by random forest using training data

%% For all observations
split = 0;          % do not split data
importSOCATv2021;   % import data
gridSOCATv2021;     % grid SOCAT data
detrendSOCATv2021;  % detrend gridded observations
loadvarsSOCATv2021; % load predictor variables
grid_all = 1;       % predict pCO2 on grid
interpolateSOCATv2021_RF_all;
                    % interpolate pCO2 by random forest using all data

%% Import other datasets
importLAND;              % Landschutzer et al. (2020)
importLAR;               % Laruelle et al. (2017)

%% Ancillary functions for calculations
pco2_uncertainty          % determine uncertainty in pCO2
calculateCO2flux          % calculate CO2 flux
mooring_omission
moorings

%% save final product as .mat and NetCDF
save_files

%% Ancillary functions for figures
Figure1;
Figure2; % need to run validation Test #2 for this to work
Figure3;
Figure4;
Figure5;
Figure6;
Figure7;
Fay_CO2_flux % calculate CO2 flux using Fay et al
Figure8;
Figure9;
FigureA1;
FigureA2;
FigureA3;
FigureA4;
FigureA5;
FigureB1;
FigureB3a;
FigureB3b;
FigureB4;

%% Differences between ECCO2 and SOCAT
% diff=...
% SOCATv2021_grid.all.salinity_ave_weighted(~isnan(SOCATv2021_grid.all.salinity_ave_weighted) & ~isnan(SOCATv2021_grid.SSS))-...
% SOCATv2021_grid.SSS(~isnan(SOCATv2021_grid.all.salinity_ave_weighted) & ~isnan(SOCATv2021_grid.SSS));
% figure; histogram(diff);
% mean(diff)
% std(diff)
% 
% diff=...
% SOCATv2021_grid.all.sst_ave_weighted(~isnan(SOCATv2021_grid.all.sst_ave_weighted) & ~isnan(SOCATv2021_grid.SST))-...
% SOCATv2021_grid.SST(~isnan(SOCATv2021_grid.all.sst_ave_weighted) & ~isnan(SOCATv2021_grid.SST));
% figure; histogram(diff);
% mean(diff)
% std(diff)

%% Percentage of filled grid cells
% Index = ~isnan(SOCATv2021_grid.pco2_RF);
% sum(~isnan(SOCATv2021_grid.all.pco2_ave_weighted(Index)))./...
% sum(~isnan(SOCATv2021_grid.SST(Index)))
