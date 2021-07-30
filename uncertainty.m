%% Uncertainty in pCO2 and CO2 Flux

% Uncertainty = u(measurement)^2 + u(grid)^2 + u(model)^2

% Create coastal and open ocean indices
Coastal_idx = ~isnan(SOCATv2021_grid.pco2_RF) & ...
    SOCATv2021_grid.distance_from_shore <= 400;
Open_idx = ~isnan(SOCATv2021_grid.pco2_RF) & ...
    SOCATv2021_grid.distance_from_shore > 400;

%% Measurement uncertainty
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
u_meas = (2*SOCATv2021.all.per_2uatm + 5*SOCATv2021.all.per_5uatm)

%% Grid uncertainty
% This is attributable to one value in a grid cell for a month not being
% representative of any location within that grid cell during any time
% within that month

% Calculate average std of all gridded values
u_grid_coast = mean(SOCATv2021_grid.all.pco2_std_unwtd(Coastal_idx),'omitnan')
u_grid_open  = mean(SOCATv2021_grid.all.pco2_std_unwtd(Open_idx),'omitnan')
%u_grid_coast = mean(SOCATv2021_grid.all.pco2_grid_uncert(Coastal_idx),'omitnan')
%u_grid_open  = mean(SOCATv2021_grid.all.pco2_grid_uncert(Open_idx),'omitnan')

%% Adjust by autocorrelation length scales
% Grid cells outside decorrelation radius
% for ln = 1:size(SOCATv2021_grid.pco2_RF_annmean,1)
%     for lt = 1:size(SOCATv2021_grid.pco2_RF_annmean,2)
%         index = ~isnan(SOCATv2021_grid.pco2_RF_annmean);
%         dist = deg2km(distance(SOCATv2021_grid.lat(lt),SOCATv2021_grid.lon(ln),...
%             SOCATv2021_grid.lat,SOCATv2021_grid.lon));
% 
%     end
% end

%% Model uncertainty
% This is attributable to the mismatch between values predicted by the
% random forest regression model and gridded values

% calculate coastal residuals and their latitude/longitude
Residuals_coast_lat = SOCATv2021_grid.latitude(Coastal_idx);
Residuals_coast_lon = SOCATv2021_grid.longitude(Coastal_idx);
Residuals_coast = SOCATv2021_grid.pco2_RF(Coastal_idx) - ...
    SOCATv2021_grid.all.pco2_ave_weighted(Coastal_idx);
Residuals_coast_lat(isnan(Residuals_coast)) = [];
Residuals_coast_lon(isnan(Residuals_coast)) = [];
Residuals_coast(isnan(Residuals_coast)) = [];

% calculate open-ocean residuals and their latitude/longitude
Residuals_open_lat = SOCATv2021_grid.latitude(Open_idx);
Residuals_open_lon = SOCATv2021_grid.longitude(Open_idx);
Residuals_open = SOCATv2021_grid.pco2_RF(Open_idx) - ...
    SOCATv2021_grid.all.pco2_ave_weighted(Open_idx);
Residuals_open_lat(isnan(Residuals_open)) = [];
Residuals_open_lon(isnan(Residuals_open)) = [];
Residuals_open(isnan(Residuals_open)) = [];

% standard deviation of coastal residuals
%u_map_coast = std(Residuals_coast)
u_map_coast = 43.36

% standard deviation of open-ocean residuals
%u_map_open = std(Residuals_open)
u_map_open = 12.41

%% Lag-1 autocorrelation with random number of points

% generate two indices for randomly selecting coastal residuals
coast_pair_index1 = randi(size(Residuals_coast,1),[1000,1]);
coast_pair_index2 = randi(size(Residuals_coast,1),[1000,1]);

% generate two indices for randomly selecting open-ocean residuals
open_pair_index1 = randi(size(Residuals_open,1),[1000,1]);
open_pair_index2 = randi(size(Residuals_open,1),[1000,1]);

% select coastal residuals and latitude/longitude
Residuals_coast_pair = ...
    [Residuals_coast(coast_pair_index1),Residuals_coast(coast_pair_index2)];
Residuals_coast_lat_pair = ...
    [Residuals_coast_lat(coast_pair_index1),Residuals_coast_lat(coast_pair_index2)];
Residuals_coast_lon_pair = ...
    [Residuals_coast_lon(coast_pair_index1),Residuals_coast_lon(coast_pair_index2)];

% select open-ocean residuals and latitude/longitude
Residuals_open_pair = ...
    [Residuals_open(open_pair_index1),Residuals_open(open_pair_index2)];
Residuals_open_lat_pair = ...
    [Residuals_open_lat(open_pair_index1),Residuals_open_lat(open_pair_index2)];
Residuals_open_lon_pair = ...
    [Residuals_open_lon(open_pair_index1),Residuals_open_lon(open_pair_index2)];

% calculate squared difference between residuals
sqr_diff_coast = ...
    (Residuals_coast_pair(:,1) - Residuals_coast_pair(:,2)).^2;
sqr_diff_open = ...
    (Residuals_open_pair(:,1) - Residuals_open_pair(:,2)).^2;

distance_coast_pair = ...
    deg2km(...
    distance(Residuals_coast_lat_pair(:,1),Residuals_coast_lon_pair(:,1),...
    Residuals_coast_lat_pair(:,2),Residuals_coast_lon_pair(:,2)));

distance_open_pair = ...
    deg2km(...
    distance(Residuals_open_lat_pair(:,1),Residuals_open_lon_pair(:,1),...
    Residuals_open_lat_pair(:,2),Residuals_open_lon_pair(:,2)));

Neff_coast = 1;
Neff_open = 1;


%% Total pCO2 uncertainty
U_tot_coast = sqrt(  u_meas.^2 + ... % measurement uncertainty
                     u_grid_coast.^2  + ... % gridding uncertainty
                     u_map_coast.^2)  % mapping uncertainty)
U_tot_open  = sqrt(  u_meas.^2 + ... % measurement uncertainty
                     u_grid_open.^2  + ... % gridding uncertainty
                     u_map_open.^2) % mapping uncertainty)

% %% Expand total uncertainty for coastal and open
% SOCATv2021_grid.pco2_RF_uncertainty = ...
%     nan(size(SOCATv2021_grid.pco2_RF));
% 
% SOCATv2021_grid.pco2_RF_uncertainty(Coastal_idx) = U_tot_coast;
% SOCATv2021_grid.pco2_RF_uncertainty(Open_idx) = U_tot_open;
% 
% %% Calculate FCO2 with added pCO2 uncertainty
% 
% SOCATv2021_grid.Fco2_RF_ERA5_minus = SOCATv2021_grid.kw_ERA5 .* ...
%     SOCATv2021_grid.k0_ERA5 .* ...
%     (SOCATv2021_grid.delpco2_RF_ERA5-SOCATv2021_grid.pco2_RF_uncertainty).*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_ERA5_none = SOCATv2021_grid.kw_ERA5 .* ...
%     SOCATv2021_grid.k0_ERA5 .* ...
%     SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_ERA5_plus = SOCATv2021_grid.kw_ERA5 .* ...
%     SOCATv2021_grid.k0_ERA5 .* ...
%     (SOCATv2021_grid.delpco2_RF_ERA5+SOCATv2021_grid.pco2_RF_uncertainty).*(1e-6); % mol C/(m^2 * month)
% 
% % Calculate standard deviation with pCO2 error
% SOCATv2021_grid.Fco2_pCO2_uncertainty = ...
% std(cat(4,SOCATv2021_grid.Fco2_RF_ERA5_minus.*12,...
%           SOCATv2021_grid.Fco2_RF_ERA5_none.*12,...
%           SOCATv2021_grid.Fco2_RF_ERA5_plus.*12),[],4);
%       
% SOCATv2021_grid.Fco2_uncertainty_pCO2_annual = ...
%     mean(SOCATv2021_grid.Fco2_pCO2_uncertainty,3,'omitnan');
% 
% %% Uncertainty in gas transfer coefficient
% 
% % Calculate Schmidt number
% SOCATv2021_grid.sch = CO2flux_Schmidt_W14(SOCATv2021_grid.SST,'CO2');
% 
% % Calculate kw using Ho et al (2006)
% SOCATv2021_grid.kw_ERA5_Ho06 = ...
%     ((0.254.*SOCATv2021_grid.wind_speed.^2.*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
%     
% % Calculate kw using Nightengale et al (2000)
% SOCATv2021_grid.kw_ERA5_Ng00 = ...
%     (((0.222.*SOCATv2021_grid.wind_speed.^2 + ...
%     0.333.*SOCATv2021_grid.wind_speed) .* ...
%     (SOCATv2021_grid.sch./600).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
% 
% % Calculate kw using Sweeney et al. (2007)
% SOCATv2021_grid.kw_ERA5_Sw07 = ...
%     ((0.27.*SOCATv2021_grid.wind_speed.^2.*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
% 
% % Calculate kw using Wanninkhof et al (2009)
% SOCATv2021_grid.kw_ERA5_Wn09 = ...
%     (((3+0.1.*SOCATv2021_grid.wind_speed+0.064.*SOCATv2021_grid.wind_speed.^2+...
%     0.011.*SOCATv2021_grid.wind_speed.^3).*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
% 
% % Calculate flux using each parameterization
% SOCATv2021_grid.Fco2_RF_ERA5_Ho06 = SOCATv2021_grid.kw_ERA5_Ho06 .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_ERA5_Ng00 = SOCATv2021_grid.kw_ERA5_Ng00 .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_ERA5_Sw07 = SOCATv2021_grid.kw_ERA5_Sw07 .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_ERA5_Wn09 = SOCATv2021_grid.kw_ERA5_Wn09 .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% 
% % Calculate standard deviation between flux calculations
% SOCATv2021_grid.Fco2_kw_uncertainty = ...
% std(cat(4,SOCATv2021_grid.Fco2_RF_ERA5_Ho06.*12,...
%           SOCATv2021_grid.Fco2_RF_ERA5_Ng00.*12,...
%           SOCATv2021_grid.Fco2_RF_ERA5_Sw07.*12),[],4);
%       
% SOCATv2021_grid.Fco2_uncertainty_kw_annual = ...
%     mean(SOCATv2021_grid.Fco2_kw_uncertainty,3,'omitnan');
% 
% %% Uncertainty in wind speed product
% 
% % Import other wind products
% 
% % Calculate Schmidt number
% SOCATv2021_grid.sch = CO2flux_Schmidt_W14(SOCATv2021_grid.SST,'CO2');
% 
% % Calculate kw using ERA5
% SOCATv2021_grid.kw_ERA5 = ...
%     ((0.254.*SOCATv2021_grid.wind_speed.^2.*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
%     
% % Calculate kw using NCEP2
% SOCATv2021_grid.kw_NCEP = ...
%     ((0.254.*SOCATv2021_grid.wind_speed_NCEP.^2.*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
% 
% % Calculate kw using CCMP
% SOCATv2021_grid.kw_CCMP = ...
%     ((0.254.*SOCATv2021_grid.wind_speed_CCMP.^2.*(SOCATv2021_grid.sch./660).^-0.5) ./ ... % cm/h
%     100).*(24.*365.2422./12); % m/month
% 
% % Calculate flux using each wind speed product
% SOCATv2021_grid.Fco2_RF_ERA5winds = SOCATv2021_grid.kw_ERA5 .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_NCEPwinds = SOCATv2021_grid.kw_NCEP .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% SOCATv2021_grid.Fco2_RF_CCMPwinds = SOCATv2021_grid.kw_CCMP .* ...
%     SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mol C/(m^2 * month)
% 
% % Calculate standard deviation between wind calculations
% SOCATv2021_grid.Fco2_wind_speed_uncertainty = ...
% std(cat(4,SOCATv2021_grid.Fco2_RF_ERA5winds.*12,...
%           SOCATv2021_grid.Fco2_RF_CCMPwinds.*12),[],4);
%       
% SOCATv2021_grid.Fco2_uncertainty_wind_annual = ...
%     mean(SOCATv2021_grid.Fco2_wind_speed_uncertainty,3,'omitnan');
% 
% %% Total FCO2 uncertainty
% 
% SOCATv2021_grid.Fco2_uncertainty_total = ...
%     sqrt(SOCATv2021_grid.Fco2_pCO2_uncertainty.^2 + ...
%          SOCATv2021_grid.Fco2_kw_uncertainty.^2 + ...
%          SOCATv2021_grid.Fco2_wind_speed_uncertainty.^2);
% 
% SOCATv2021_grid.Fco2_uncertainty_total_annual = ...
%     mean(SOCATv2021_grid.Fco2_uncertainty_total,3,'omitnan');
% 
% % Weighted mean over domain
% area_weights_coastal = SOCATv2021_grid.area_km2;
% area_weights_open = SOCATv2021_grid.area_km2;
% area_weights_coastal(isnan(SOCATv2021_grid.pco2_RF(:,:,1))) = NaN;
% area_weights_coastal(~Coastal_idx(:,:,1)) = NaN;
% area_weights_open(isnan(SOCATv2021_grid.pco2_RF(:,:,1))) = NaN;
% area_weights_open(~Open_idx(:,:,1)) = NaN;
% 
% Fco2_uncertainty_pCO2_coastal =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_pCO2_annual.*area_weights_coastal,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_coastal,1,'omitnan'),2,'omitnan'))...
%     )
% Fco2_uncertainty_pCO2_open =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_pCO2_annual.*area_weights_open,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_open,1,'omitnan'),2,'omitnan'))...
%     )
% 
% Fco2_uncertainty_kw_coastal =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_kw_annual.*area_weights_coastal,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_coastal,1,'omitnan'),2,'omitnan'))...
%     )
% Fco2_uncertainty_kw_open =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_kw_annual.*area_weights_open,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_open,1,'omitnan'),2,'omitnan'))...
%     )
% 
% Fco2_uncertainty_wind_coastal =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_wind_annual.*area_weights_coastal,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_coastal,1,'omitnan'),2,'omitnan'))...
%     )
% Fco2_uncertainty_wind_open =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_wind_annual.*area_weights_open,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_open,1,'omitnan'),2,'omitnan'))...
%     )
% 
% Fco2_uncertainty_tot_coastal =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_total_annual.*area_weights_coastal,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_coastal,1,'omitnan'),2,'omitnan'))...
%     )
% Fco2_uncertainty_tot_open =...
%     mean(...
%     squeeze(sum(sum(SOCATv2021_grid.Fco2_uncertainty_total_annual.*area_weights_open,1,'omitnan'),2,'omitnan'))./...
%     squeeze(sum(sum(area_weights_open,1,'omitnan'),2,'omitnan'))...
%     )
% 
% % pcolor(SOCATv2021_grid.lon,SOCATv2021_grid.lat,SOCATv2021_grid.Fco2_uncertainty_total(:,:,122))