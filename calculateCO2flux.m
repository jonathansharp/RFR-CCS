%% Calculate CO2 Flux for RF pCO2 product

%% Calculate flux with Humphreys toolbox function

[SOCATv2020_grid.Fco2_RF_ERA5_hump, ~, ~, ...
    SOCATv2020_grid.delpco2_RF_ERA5_hump, ...
    SOCATv2020_grid.kw_ERA5_hump, ... % cm/hr
    SOCATv2020_grid.k0_ERA5_hump] = ... % mol/(l*atm)
CO2flux(SOCATv2020_grid.SST,SOCATv2020_grid.SSS,...
    SOCATv2020_grid.wind_speed,SOCATv2020_grid.pCO2_atm, ...
    0,SOCATv2020_grid.pco2_RF,0,'w14'); % umol/(m^2*hr)

SOCATv2020_grid.Fco2_RF_ERA5_hump = ...
   SOCATv2020_grid.Fco2_RF_ERA5_hump.*(24.*365.25).*(1e-6); % mol/(m^2*yr)

%% Calculate flux manually

% Calculate gas transfer velocity
SOCATv2020_grid.kw_ERA5 = kgas(SOCATv2020_grid.wind_speed,660,'W14'); % m/s

% Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
TempK100  = (SOCATv2020_grid.SST+273.15)./100;
lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + SOCATv2020_grid.SSS .*...
    (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
SOCATv2020_grid.k0_ERA5 = exp(lnK0); % mol/kg-SW/atm

% Calculate density
SOCATv2020_grid.SA      = gsw_SA_from_SP(SOCATv2020_grid.SSS,5,SOCATv2020_grid.longitude,SOCATv2020_grid.latitude); % Calculate absolute salinity
SOCATv2020_grid.CT      = gsw_CT_from_t(SOCATv2020_grid.SA,SOCATv2020_grid.SST,5); % Calculate conservative temp
SOCATv2020_grid.DENSITY = gsw_rho(SOCATv2020_grid.SA,SOCATv2020_grid.CT,5); % Calculate in situ density

% Calculate delta pCO2
SOCATv2020_grid.delpco2_RF_ERA5 = SOCATv2020_grid.pco2_RF - SOCATv2020_grid.pCO2_atm; % uatm

% Calculate flux
SOCATv2020_grid.Fco2_RF_ERA5 = SOCATv2020_grid.kw_ERA5 .* SOCATv2020_grid.k0_ERA5 .* (SOCATv2020_grid.delpco2_RF_ERA5.*1e-6); % m/s
SOCATv2020_grid.Fco2_RF_ERA5 = SOCATv2020_grid.Fco2_RF_ERA5.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
SOCATv2020_grid.Fco2_RF_ERA5 = SOCATv2020_grid.Fco2_RF_ERA5.*SOCATv2020_grid.DENSITY; % mol/(m^2*yr)

%% Calculate flux manually with U^2 held constant and with delta pCO2 held constant

% Calculate gas transfer velocity (U^2 held constant)
Uconst = sqrt(mean(SOCATv2020_grid.wind_speed.^2,3,'omitnan'));
SOCATv2020_grid.kw_ERA5_Uconst = kgas(Uconst,660,'W14'); % m/s
SOCATv2020_grid.kw_ERA5_Uconst = repmat(SOCATv2020_grid.kw_ERA5_Uconst,1,1,size(SOCATv2020_grid.month_since_1998,1));

% Calculate delta pCO2 (delta pCO2 held constant)
SOCATv2020_grid.delpco2_RF_ERA5_PCO2const = mean(SOCATv2020_grid.pco2_RF - SOCATv2020_grid.pCO2_atm,3,'omitnan'); % uatm

% Calculate flux (U^2 held constant)
SOCATv2020_grid.Fco2_RF_ERA5_Uconst = SOCATv2020_grid.kw_ERA5_Uconst .* SOCATv2020_grid.k0_ERA5 .* (SOCATv2020_grid.delpco2_RF_ERA5.*1e-6); % m/s
SOCATv2020_grid.Fco2_RF_ERA5_Uconst = SOCATv2020_grid.Fco2_RF_ERA5_Uconst.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
SOCATv2020_grid.Fco2_RF_ERA5_Uconst = SOCATv2020_grid.Fco2_RF_ERA5_Uconst.*SOCATv2020_grid.DENSITY; % mol/(m^2*yr)

% Calculate flux (delta pCO2 held constant)
SOCATv2020_grid.Fco2_RF_ERA5_PCO2const = SOCATv2020_grid.kw_ERA5 .* SOCATv2020_grid.k0_ERA5 .* (SOCATv2020_grid.delpco2_RF_ERA5_PCO2const.*1e-6); % m/s
SOCATv2020_grid.Fco2_RF_ERA5_PCO2const = SOCATv2020_grid.Fco2_RF_ERA5_PCO2const.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
SOCATv2020_grid.Fco2_RF_ERA5_PCO2const = SOCATv2020_grid.Fco2_RF_ERA5_PCO2const.*SOCATv2020_grid.DENSITY; % mol/(m^2*yr)

%% Determine RMSEs (averaged over one annual cycle) associated with holding U^2 constant and holding delta pCO2 constant

% Average CO2 flux components over one annual cycle
for m=1:12
    SOCATv2020_grid.Fco2_RF_ERA5_avg(:,:,m) = mean(SOCATv2020_grid.Fco2_RF_ERA5(:,:,m:12:end),3,'omitnan');
    SOCATv2020_grid.Fco2_RF_ERA5_Uconst_avg(:,:,m) = mean(SOCATv2020_grid.Fco2_RF_ERA5_Uconst(:,:,m:12:end),3,'omitnan');
    SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_avg(:,:,m) = mean(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const(:,:,m:12:end),3,'omitnan');
end

% Determine RMSEs for each component
SOCATv2020_grid.RMSE_RF_ERA5_Uconst_avg = ...
    sqrt(mean((SOCATv2020_grid.Fco2_RF_ERA5_Uconst_avg - SOCATv2020_grid.Fco2_RF_ERA5_avg).^2,3,'omitnan'));
SOCATv2020_grid.RMSE_RF_ERA5_PCO2const_avg = ...
    sqrt(mean((SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_avg - SOCATv2020_grid.Fco2_RF_ERA5_avg).^2,3,'omitnan'));

% Determine RMSE ratio
SOCATv2020_grid.RMSE_ratio_RF_ERA5_avg = SOCATv2020_grid.RMSE_RF_ERA5_Uconst_avg./SOCATv2020_grid.RMSE_RF_ERA5_PCO2const_avg;

%% Plot flux for a given grid cell over time
idxlon = 23; % 225.625
idxlat = 81; % 35.125
figure; hold on;
plot(1:12,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_avg(idxlon,idxlat,:)),'k');
plot(1:12,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_Uconst_avg(idxlon,idxlat,:)),'r');
plot(1:12,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_avg(idxlon,idxlat,:)),'b');

figure; scatter(squeeze(SOCATv2020_grid.Fco2_RF_ERA5_avg(idxlon,idxlat,:)),squeeze(SOCATv2020_grid.Fco2_RF_ERA5_Uconst_avg(idxlon,idxlat,:)));
figure; scatter(squeeze(SOCATv2020_grid.Fco2_RF_ERA5_avg(idxlon,idxlat,:)),squeeze(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_avg(idxlon,idxlat,:)));

%% Determine RMSEs (from detrended time series) associated with holding U^2 constant and holding delta pCO2 constant

% Detrend CO2 flux components
SOCATv2020_grid.Fco2_RF_ERA5_dtr = nan(size(SOCATv2020_grid.pco2_RF));
SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr = nan(size(SOCATv2020_grid.pco2_RF));
SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr = nan(size(SOCATv2020_grid.pco2_RF));
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2020_grid.Fco2_RF_ERA5(a,b,:)))) == 0
            [~,yr] = leastsq2(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5(a,b,:)),0,0,0);
            SOCATv2020_grid.Fco2_RF_ERA5_dtr(a,b,:) = yr;
            [~,yr] = leastsq2(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_Uconst(a,b,:)),0,0,0);
            SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr(a,b,:) = yr;
            [~,yr] = leastsq2(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const(a,b,:)),0,0,0);
            SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr(a,b,:) = yr;
        else
            SOCATv2020_grid.Fco2_RF_ERA5_dtr(a,b,:) = NaN;
            SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr(a,b,:) = NaN;
            SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr(a,b,:) = NaN;
        end
    end
end

% Determine RMSEs for each component
SOCATv2020_grid.RMSE_RF_ERA5_Uconst_dtr = ...
    sqrt(mean((SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr - SOCATv2020_grid.Fco2_RF_ERA5_dtr).^2,3,'omitnan'));
SOCATv2020_grid.RMSE_RF_ERA5_PCO2const_dtr = ...
    sqrt(mean((SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr - SOCATv2020_grid.Fco2_RF_ERA5_dtr).^2,3,'omitnan'));

% Determine RMSE ratio
SOCATv2020_grid.RMSE_ratio_RF_ERA5_dtr = SOCATv2020_grid.RMSE_RF_ERA5_Uconst_dtr./SOCATv2020_grid.RMSE_RF_ERA5_PCO2const_dtr;

%% Determine correlation between wind speed and delta pCO2

% Detrend components
SOCATv2020_grid.wind_speed_dtr = nan(size(SOCATv2020_grid.wind_speed));
SOCATv2020_grid.delpco2_RF_ERA5_dtr = nan(size(SOCATv2020_grid.delpco2_RF_ERA5));
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2020_grid.delpco2_RF_ERA5(a,b,:)))) == 0
            [~,yr] = leastsq2(1:264,squeeze(SOCATv2020_grid.wind_speed(a,b,:)),0,0,0);
            SOCATv2020_grid.wind_speed_dtr(a,b,:) = yr;
            [~,yr] = leastsq2(1:264,squeeze(SOCATv2020_grid.delpco2_RF_ERA5(a,b,:)),0,0,0);
            SOCATv2020_grid.delpco2_RF_ERA5_dtr(a,b,:) = yr;
        else
            SOCATv2020_grid.wind_speed_dtr(a,b,:) = NaN;
            SOCATv2020_grid.delpco2_RF_ERA5_dtr(a,b,:) = NaN;
        end
    end
end

% Determine correlation between wind speed and delta pCO2
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2020_grid.delpco2_RF_ERA5(a,b,:)))) == 0
            SOCATv2020_grid.U_vs_delpCO2_corr(a,b,:) = ...
                corr(squeeze(SOCATv2020_grid.wind_speed_dtr(a,b,:)),...
                     squeeze(SOCATv2020_grid.delpco2_RF_ERA5_dtr(a,b,:)));
        else
            SOCATv2020_grid.U_vs_delpCO2_corr(a,b,:) = NaN;
        end
    end
end

%% Determine season of max wind speed, pCO2, and temp
% Determine seasonal means
SOCATv2020_grid.delpco2_RF_ERA5_seasmean(:,:,1) = ...
    mean(cat(3,SOCATv2020_grid.delpco2_RF_ERA5(:,:,1:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,2:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,12:12:end)),3,'omitnan');
SOCATv2020_grid.delpco2_RF_ERA5_seasmean(:,:,2) = ...
    mean(cat(3,SOCATv2020_grid.delpco2_RF_ERA5(:,:,3:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,4:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,5:12:end)),3,'omitnan');
SOCATv2020_grid.delpco2_RF_ERA5_seasmean(:,:,3) = ...
    mean(cat(3,SOCATv2020_grid.delpco2_RF_ERA5(:,:,6:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,7:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,8:12:end)),3,'omitnan');
SOCATv2020_grid.delpco2_RF_ERA5_seasmean(:,:,4) = ...
    mean(cat(3,SOCATv2020_grid.delpco2_RF_ERA5(:,:,9:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,10:12:end),...
               SOCATv2020_grid.delpco2_RF_ERA5(:,:,11:12:end)),3,'omitnan');
SOCATv2020_grid.wind_speed_seasmean(:,:,1) = ...
    mean(cat(3,SOCATv2020_grid.wind_speed(:,:,1:12:end),...
               SOCATv2020_grid.wind_speed(:,:,2:12:end),...
               SOCATv2020_grid.wind_speed(:,:,12:12:end)),3,'omitnan');
SOCATv2020_grid.wind_speed_seasmean(:,:,2) = ...
    mean(cat(3,SOCATv2020_grid.wind_speed(:,:,3:12:end),...
               SOCATv2020_grid.wind_speed(:,:,4:12:end),...
               SOCATv2020_grid.wind_speed(:,:,5:12:end)),3,'omitnan');
SOCATv2020_grid.wind_speed_seasmean(:,:,3) = ...
    mean(cat(3,SOCATv2020_grid.wind_speed(:,:,6:12:end),...
               SOCATv2020_grid.wind_speed(:,:,7:12:end),...
               SOCATv2020_grid.wind_speed(:,:,8:12:end)),3,'omitnan');
SOCATv2020_grid.wind_speed_seasmean(:,:,4) = ...
    mean(cat(3,SOCATv2020_grid.wind_speed(:,:,9:12:end),...
               SOCATv2020_grid.wind_speed(:,:,10:12:end),...
               SOCATv2020_grid.wind_speed(:,:,11:12:end)),3,'omitnan');
SOCATv2020_grid.SST_seasmean(:,:,1) = ...
    mean(cat(3,SOCATv2020_grid.SST(:,:,1:12:end),...
               SOCATv2020_grid.SST(:,:,2:12:end),...
               SOCATv2020_grid.SST(:,:,12:12:end)),3,'omitnan');
SOCATv2020_grid.SST_seasmean(:,:,2) = ...
    mean(cat(3,SOCATv2020_grid.SST(:,:,3:12:end),...
               SOCATv2020_grid.SST(:,:,4:12:end),...
               SOCATv2020_grid.SST(:,:,5:12:end)),3,'omitnan');
SOCATv2020_grid.SST_seasmean(:,:,3) = ...
    mean(cat(3,SOCATv2020_grid.SST(:,:,6:12:end),...
               SOCATv2020_grid.SST(:,:,7:12:end),...
               SOCATv2020_grid.SST(:,:,8:12:end)),3,'omitnan');
SOCATv2020_grid.SST_seasmean(:,:,4) = ...
    mean(cat(3,SOCATv2020_grid.SST(:,:,9:12:end),...
               SOCATv2020_grid.SST(:,:,10:12:end),...
               SOCATv2020_grid.SST(:,:,11:12:end)),3,'omitnan');
SOCATv2020_grid.CHL_seasmean(:,:,1) = ...
    mean(cat(3,SOCATv2020_grid.CHL(:,:,1:12:end),...
               SOCATv2020_grid.CHL(:,:,2:12:end),...
               SOCATv2020_grid.CHL(:,:,12:12:end)),3,'omitnan');
SOCATv2020_grid.CHL_seasmean(:,:,2) = ...
    mean(cat(3,SOCATv2020_grid.CHL(:,:,3:12:end),...
               SOCATv2020_grid.CHL(:,:,4:12:end),...
               SOCATv2020_grid.CHL(:,:,5:12:end)),3,'omitnan');
SOCATv2020_grid.CHL_seasmean(:,:,3) = ...
    mean(cat(3,SOCATv2020_grid.CHL(:,:,6:12:end),...
               SOCATv2020_grid.CHL(:,:,7:12:end),...
               SOCATv2020_grid.CHL(:,:,8:12:end)),3,'omitnan');
SOCATv2020_grid.CHL_seasmean(:,:,4) = ...
    mean(cat(3,SOCATv2020_grid.CHL(:,:,9:12:end),...
               SOCATv2020_grid.CHL(:,:,10:12:end),...
               SOCATv2020_grid.CHL(:,:,11:12:end)),3,'omitnan');
           
SOCATv2020_grid.seasmax_delpCO2 = nan(size(SOCATv2020_grid.delpco2_RF_ERA5,1),size(SOCATv2020_grid.delpco2_RF_ERA5,2));
SOCATv2020_grid.seasmax_wind_speed = nan(size(SOCATv2020_grid.wind_speed,1),size(SOCATv2020_grid.wind_speed,2));
SOCATv2020_grid.seasmax_SST = nan(size(SOCATv2020_grid.SST,1),size(SOCATv2020_grid.SST,2));
SOCATv2020_grid.seasmax_CHL = nan(size(SOCATv2020_grid.CHL,1),size(SOCATv2020_grid.CHL,2));
for a = 1:size(SOCATv2020_grid.delpco2_RF_ERA5_seasmean,1)
    for b = 1:size(SOCATv2020_grid.delpco2_RF_ERA5_seasmean,2)
        if sum(isnan(squeeze(SOCATv2020_grid.delpco2_RF_ERA5_seasmean(a,b,:)))) == 4
            SOCATv2020_grid.seasmax_delpCO2(a,b) = NaN;
            SOCATv2020_grid.seasmax_wind_speed(a,b) = NaN;
            SOCATv2020_grid.seasmax_SST(a,b) = NaN;
            SOCATv2020_grid.seasmax_CHL(a,b) = NaN;
        else
            SOCATv2020_grid.seasmax_delpCO2(a,b) = find(abs(squeeze(SOCATv2020_grid.delpco2_RF_ERA5_seasmean(a,b,:))) == max(abs(squeeze(SOCATv2020_grid.delpco2_RF_ERA5_seasmean(a,b,:)))));
            SOCATv2020_grid.seasmax_wind_speed(a,b) = find(squeeze(SOCATv2020_grid.wind_speed_seasmean(a,b,:)) == max(squeeze(SOCATv2020_grid.wind_speed_seasmean(a,b,:))));
            SOCATv2020_grid.seasmax_SST(a,b) = find(squeeze(SOCATv2020_grid.SST_seasmean(a,b,:)) == max(squeeze(SOCATv2020_grid.SST_seasmean(a,b,:))));
            SOCATv2020_grid.seasmax_CHL(a,b) = find(squeeze(SOCATv2020_grid.CHL_seasmean(a,b,:)) == max(squeeze(SOCATv2020_grid.CHL_seasmean(a,b,:))));
        end
    end        
end

%% Plot flux for a given grid cell over time
idxlon = 23; % 225.625
idxlat = 81; % 35.125
figure; hold on;
plot(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_dtr(idxlon,idxlat,:)),'k');
plot(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr(idxlon,idxlat,:)),'r');
plot(1:264,squeeze(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr(idxlon,idxlat,:)),'b');

figure; scatter(squeeze(SOCATv2020_grid.Fco2_RF_ERA5_dtr(idxlon,idxlat,:)),squeeze(SOCATv2020_grid.Fco2_RF_ERA5_Uconst_dtr(idxlon,idxlat,:)));
figure; scatter(squeeze(SOCATv2020_grid.Fco2_RF_ERA5_dtr(idxlon,idxlat,:)),squeeze(SOCATv2020_grid.Fco2_RF_ERA5_PCO2const_dtr(idxlon,idxlat,:)));

%% Calculate CO2 Flux for Landschutzer pCO2 product
% %% Calculate flux manually
% % Gas transfer velocity
% for m = 1:12
%     LAND.wind_speed(:,:,m) = mean(SOCATv2020_grid.wind_speed(:,:,m:12:216),3,'omitnan');
%     LAND.SST(:,:,m) = mean(SOCATv2020_grid.SST(:,:,m:12:216),3,'omitnan');
%     LAND.SSS(:,:,m) = mean(SOCATv2020_grid.SSS(:,:,m:12:216),3,'omitnan');
%     LAND.pCO2_atm(:,:,m) = mean(SOCATv2020_grid.pCO2_atm(:,:,m:12:216),3,'omitnan');
% end
% LAND.kw_ERA5 = kgas(LAND.wind_speed,660,'W14'); % m/s
% 
% % Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
% TempK100  = (LAND.SST+273.15)./100;
% lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + LAND.SSS .*...
%     (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
% LAND.k0_ERA5 = exp(lnK0); % mol/kg-SW/atm
% 
% % Calculate density
% LAND.SA      = gsw_SA_from_SP(LAND.SSS,5,LAND.longitude,LAND.latitude); % Calculate absolute salinity
% LAND.CT      = gsw_CT_from_t(LAND.SA,LAND.SST,5); % Calculate conservative temp
% LAND.DENSITY = gsw_rho(LAND.SA,LAND.CT,5); % Calculate in situ density
% 
% % Calculate delta pCO2
% LAND.delpco2_ERA5 = LAND.pCO2 - LAND.pCO2_atm; % uatm
% 
% % Calculate flux
% 
% LAND.Fco2_ERA5 = LAND.kw_ERA5 .* LAND.k0_ERA5 .* (LAND.delpco2_ERA5.*1e-6); % m/s
% LAND.Fco2_ERA5 = LAND.Fco2_ERA5.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
% LAND.Fco2_ERA5 = LAND.Fco2_ERA5.*LAND.DENSITY; % mol/(m^2*yr)

