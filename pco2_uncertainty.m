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
n_meas = sum(~isnan(SOCATv2021.all.pCO2))

%% Grid uncertainty
% This is attributable to one value in a grid cell for a month not being
% representative of any location within that grid cell during any time
% within that month

% Calculate average std of all gridded values
u_grid_coast = mean(SOCATv2021_grid.all.pco2_grid_uncert(Coastal_idx),'omitnan')
n_grid_coast = sum(~isnan(SOCATv2021_grid.all.pco2_grid_uncert(Coastal_idx)))
u_grid_open  = mean(SOCATv2021_grid.all.pco2_grid_uncert(Open_idx),'omitnan')
n_grid_open = sum(~isnan(SOCATv2021_grid.all.pco2_grid_uncert(Open_idx)))

% % Calculate average number of points outside of decorrelation length scale
% % for open-ocean (400km) observations
% idx_open_temp = ~isnan(SOCATv2021_grid.pco2_RF(:,:,1)) & ...
%     SOCATv2021_grid.distance_from_shore(:,:,1) > 400;
% idx_open_temp = idx_open_temp(:);
% lat_temp = SOCATv2021_grid.lat(:); lat_temp = lat_temp(idx_open_temp);
% lon_temp = SOCATv2021_grid.lon(:); lon_temp = lon_temp(idx_open_temp);
% Neff_open = nan(size(lat_temp));
% for n = 1:length(lat_temp)
%     dist = distance(lat_temp,lon_temp,lat_temp(n),lon_temp(n));
%     dist = deg2km(dist);
%     Neff_open(n) = sum(sum(dist > 400));
% end
% Neff_open_mean = mean(Neff_open);
% 
% % Calculate average number of points outside of decorrelation length scale
% % for coastal (50km) observations
% idx_coast_temp = ~isnan(SOCATv2021_grid.pco2_RF(:,:,1)) & ...
%     SOCATv2021_grid.distance_from_shore(:,:,1) > 50;
% idx_coast_temp = idx_coast_temp(:);
% lat_temp = SOCATv2021_grid.lat(:); lat_temp = lat_temp(idx_coast_temp);
% lon_temp = SOCATv2021_grid.lon(:); lon_temp = lon_temp(idx_coast_temp);
% Neff_coast = nan(size(lat_temp));
% for n = 1:length(lat_temp)
%     dist = distance(lat_temp,lon_temp,lat_temp(n),lon_temp(n));
%     dist = deg2km(dist);
%     Neff_coast(n) = sum(sum(dist > 400));
% end
% Neff_coast_mean = mean(Neff_coast);

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


