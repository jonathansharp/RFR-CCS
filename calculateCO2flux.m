%% Calculate CO2 Flux for RF pCO2 product

% Calculate gas transfer velocity with kgas function
SOCATv2021_grid.sch = CO2flux_Schmidt_W14(SOCATv2021_grid.SST,'CO2');
SOCATv2021_grid.kw_ERA5 = kgas(SOCATv2021_grid.wind_speed,SOCATv2021_grid.sch,'Ho06'); % m/s
SOCATv2021_grid.kw_ERA5 = SOCATv2021_grid.kw_ERA5.*(60.*60.*24.*365.25./12); % m/month
SOCATv2021_grid.kw_ERA5_hourly = SOCATv2021_grid.kw_ERA5./(24.*365.25./12); % m/hr

% Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
TempK100  = (SOCATv2021_grid.SST+273.15)./100;
lnK0 = -58.0931 + 90.5069 ./ TempK100 + 22.2940 .* log(TempK100) + SOCATv2021_grid.SSS .*...
    (0.027766 - 0.025888 .* TempK100 + 0.0050578 .* TempK100 .^2);
SOCATv2021_grid.k0_ERA5 = exp(lnK0).*1e3; % mmol/(l*atm)
SOCATv2021_grid.k0_ERA5 = SOCATv2021_grid.k0_ERA5.*1000; % mmol/(m^3 * atm)

% Calculate delta pCO2
SOCATv2021_grid.delpco2_RF_ERA5 = SOCATv2021_grid.pco2_RF - SOCATv2021_grid.pCO2_atm; % uatm

% Calculate flux

SOCATv2021_grid.Fco2_RF_ERA5 = SOCATv2021_grid.kw_ERA5 .* ...
    SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mmol C/(m^2 * month)

SOCATv2021_grid.Fco2_RF_ERA5_hourly = SOCATv2021_grid.kw_ERA5_hourly .* ...
    SOCATv2021_grid.k0_ERA5 .* SOCATv2021_grid.delpco2_RF_ERA5.*(1e-6); % mmol C/(m^2 * hr)

%% Calculate averages
% Average CO2 flux over one annual cycle
for m=1:12
    SOCATv2021_grid.Fco2_RF_ERA5_clim(:,:,m) = mean(SOCATv2021_grid.Fco2_RF_ERA5(:,:,m:12:end),3,'omitnan');
end
% Average CO2 flux for each year
for y=1:22
    SOCATv2021_grid.Fco2_RF_ERA5_annual(:,:,y) = mean(SOCATv2021_grid.Fco2_RF_ERA5(:,:,(y-1)*12+1:(y-1)*12+12),3,'omitnan');
end
% Annual mean over time period
SOCATv2021_grid.Fco2_RF_ERA5_mean = mean(SOCATv2021_grid.Fco2_RF_ERA5_clim,3,'omitnan');

%% Calculate standard deviations
% Monthly CO2 flux standard deviation
for m=1:12
    SOCATv2021_grid.Fco2_monthly_std(:,:,m) = ...
        std(SOCATv2021_grid.Fco2_RF_ERA5(:,:,m:12:end),[],3);
end
% Yearly CO2 flux standard deviation
for y=1:22
    SOCATv2021_grid.Fco2_yearly_std(:,:,y) = ...
        std(SOCATv2021_grid.Fco2_RF_ERA5(:,:,(y-1)*12+1:(y-1)*12+12),[],3);
end

%% Calculate amplitude
SOCATv2021_grid.Fco2_RF_ERA5_amp = ...
    max(SOCATv2021_grid.Fco2_RF_ERA5_clim,[],3)-...
    min(SOCATv2021_grid.Fco2_RF_ERA5_clim,[],3);

%% Calculate domain means
area_weights = SOCATv2021_grid.area_km2;
area_weights = repmat(area_weights,1,1,size(SOCATv2021_grid.Fco2_RF_ERA5,3));
area_weights(isnan(SOCATv2021_grid.Fco2_RF_ERA5)) = NaN;
SOCATv2021_grid.Fco2_RF_dom_mean = ...
    squeeze(sum(sum(SOCATv2021_grid.Fco2_RF_ERA5.*area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
SOCATv2021_grid.Fco2_RF_dom_std = nan(size(SOCATv2021_grid.month_since_1998,1),1);
for a = 1:size(SOCATv2021_grid.month_since_1998,1)
    FCO2temp = SOCATv2021_grid.Fco2_RF_ERA5(:,:,a);
    weighttemp = area_weights(:,:,a);
    SOCATv2021_grid.Fco2_RF_dom_std(a) = std(FCO2temp(:),weighttemp(:),'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate flux manually (Laruelle)

% Extract gas transfer velocity from 1998-2015
LAR.kw_ERA5 = SOCATv2021_grid.kw_ERA5(:,:,1:216);
LAR.kw_ERA5_hourly = SOCATv2021_grid.kw_ERA5_hourly(:,:,1:216);

% Extract K0 from 1998-2015
LAR.k0_ERA5 = SOCATv2021_grid.k0_ERA5(:,:,1:216);

% Calculate delta pCO2 from 1998-2015
LAR.delpco2_ERA5 = LAR.pCO2 - SOCATv2021_grid.pCO2_atm(:,:,1:216); % uatm

% Calculate flux
LAR.Fco2_ERA5 = LAR.kw_ERA5 .* LAR.k0_ERA5 .* LAR.delpco2_ERA5.*(1e-6); % mol C/(m^2 * month)
LAR.Fco2_ERA5_hourly = LAR.kw_ERA5_hourly .* LAR.k0_ERA5 .* LAR.delpco2_ERA5.*(1e-6); % mol C/(m^2 * hr)

%% Calculate averages
% Average CO2 flux over one annual cycle
for m=1:12
    LAR.Fco2_ERA5_clim(:,:,m) = mean(LAR.Fco2_ERA5(:,:,m:12:end),3,'omitnan');
end
% Average CO2 flux for each year
for y=1:18
    LAR.Fco2_ERA5_annual(:,:,y) = mean(LAR.Fco2_ERA5(:,:,(y-1)*12+1:(y-1)*12+12),3,'omitnan');
end

%% Calculate standard deviations
% Monthly CO2 flux standard deviation
for m=1:12
    LAR.Fco2_monthly_std(:,:,m) = ...
        std(LAR.Fco2_ERA5(:,:,m:12:end),[],3);
end
% Yearly CO2 flux standard deviation
for y=1:18
    LAR.Fco2_yearly_std(:,:,y) = ...
        std(LAR.Fco2_ERA5(:,:,(y-1)*12+1:(y-1)*12+12),[],3);
end

%% Calculate amplitude
LAR.Fco2_ERA5_amp = ...
    max(LAR.Fco2_ERA5_clim,[],3)-...
    min(LAR.Fco2_ERA5_clim,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate flux manually (Landschutzer)

% Calculate average gas transfer velocity from 1998-2015
LAND.kw_ERA5 = nan(size(LAND.pCO2));
for m=1:12
LAND.kw_ERA5(:,:,m) = mean(SOCATv2021_grid.kw_ERA5(:,:,m:12:18*12),3,'omitnan'); % m/month
end

% Calculate average K0 from 1998-2015
LAND.k0_ERA5 = nan(size(LAND.pCO2));
for m=1:12
LAND.k0_ERA5(:,:,m) = mean(SOCATv2021_grid.k0_ERA5(:,:,m:12:18*12),3,'omitnan'); % mol/(m^3 * atm)
end

% Calculate average atmospheric pCO2 from 1998-2015
LAND.pCO2_atm = nan(size(LAND.pCO2));
for m=1:12
LAND.pCO2_atm(:,:,m) = mean(SOCATv2021_grid.pCO2_atm(:,:,m:12:18*12),3,'omitnan'); % uatm
end

% Calculate delta pCO2
LAND.delpco2_ERA5 = LAND.pCO2 - LAND.pCO2_atm; % uatm

% Calculate flux
LAND.Fco2_ERA5 = LAND.kw_ERA5 .* LAND.k0_ERA5 .* LAND.delpco2_ERA5.*(1e-6); % mol C/(m^2 * month)

%% Calculate amplitude
LAND.Fco2_ERA5_amp = ...
    max(LAND.Fco2_ERA5,[],3)-...
    min(LAND.Fco2_ERA5,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
