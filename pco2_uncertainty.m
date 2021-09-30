%% Uncertainty in pCO2

% Create coastal and open ocean indices
Coastal_idx = ~isnan(SOCATv2021_grid.pco2_RF) & ...
    SOCATv2021_grid.distance_from_shore <= 400;
Open_idx = ~isnan(SOCATv2021_grid.pco2_RF) & ...
    SOCATv2021_grid.distance_from_shore > 400;

%% Measurement uncertainty in pCO2
% Measurement uncertainty is either 2 uatm or 5 uatm

% Convert dataset flags to numeric:
SOCATv2021.all.flagn = SOCATv2021.all.flag;
SOCATv2021.all.flagn(strcmp(SOCATv2021.all.flag,'A')) = {1};
SOCATv2021.all.flagn(strcmp(SOCATv2021.all.flag,'B')) = {2};
SOCATv2021.all.flagn(strcmp(SOCATv2021.all.flag,'C')) = {3};
SOCATv2021.all.flagn(strcmp(SOCATv2021.all.flag,'D')) = {4};
SOCATv2021.all.flagn = cell2mat(SOCATv2021.all.flagn);

% Percent with estimated uncertainty of 2 uatm vs. 5 uatm:
SOCATv2021.all.per_2uatm = ...
    sum(SOCATv2021.all.flagn == 1 | SOCATv2021.all.flagn == 2)./ ...
    size(SOCATv2021.all.flagn,1);
SOCATv2021.all.per_5uatm = 1 - SOCATv2021.all.per_2uatm;

% Average uncertainty:
u_meas = 2*SOCATv2021.all.per_2uatm + 5*SOCATv2021.all.per_5uatm

%% Grid uncertainty
% This is attributable to one value in a grid cell for a month not being
% representative of any location within that grid cell during any time
% within that month

% Calculate average std of all gridded values
u_grid_coast = mean(SOCATv2021_grid.all.pco2_grid_uncert(Coastal_idx),'omitnan')
u_grid_open  = mean(SOCATv2021_grid.all.pco2_grid_uncert(Open_idx),'omitnan')

%% Map uncertainty in pCO2
% This is attributable to the mismatch metween values predicted by the
% random forest regression model and gridded values

% Compare predicted values to test values (coastal)
% u_map_coast = std(SOCATv2021_grid.pco2_RF(Coastal_idx) - ...
%     SOCATv2021_grid.all.pco2_ave_weighted(Coastal_idx),[],'omitnan')

% Compare predicted values to test values (open)
% u_map_open = std(SOCATv2021_grid.pco2_RF(~Coastal_idx) - ...
%     SOCATv2021_grid.all.pco2_ave_weighted(~Coastal_idx),[],'omitnan')

u_map_coast = 35.27 % determined as average RMSE from independent tests #1 and #2 for coastal ocean 
u_map_open = 4.39 % determined as average RMSE from independent tests #1 and #2 for open ocean

%% Combined uncertainty in pCO2

U_tot_coast = sqrt(u_meas.^2 + u_grid_coast.^2 + u_map_coast.^2)
U_tot_open  = sqrt(u_meas.^2 + u_grid_open.^2 + u_map_open.^2)


