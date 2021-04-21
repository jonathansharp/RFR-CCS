%% SOCATv2020_grid
% This script creates a climatology of surface pCO2 for eastern Pacific
% region that includes the California Current Large Marine Ecosystem

%disp('Loading my gridded SOCAT file truncated to CCS region');
%load('/Volumes/2TB Hard Drive/SOCAT/socatv2020_jongrid.mat');

%% Determine months and years for gridded products
disp('Converting months since 1998 to years and months of year');
date = datetime([repmat(1998,size(SOCATv2020_grid.month_since_1998,1),1),...
              SOCATv2020_grid.month_since_1998,...
              ones(size(SOCATv2020_grid.month_since_1998,1),1)]);
date = datevec(date);
SOCATv2020_grid.year = date(:,1);
SOCATv2020_grid.month = date(:,2);
clear date

%% Extend SOCAT lat, lon, and month through time in SOCATv2020_combined structure
SOCATv2020_grid.month_since_1998_3D = ...
    repmat(permute(SOCATv2020_grid.month_since_1998,[3 2 1]),...
    size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),1);
SOCATv2020_grid.year = ...
    repmat(permute(SOCATv2020_grid.year,[3 2 1]),...
    size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),1);
SOCATv2020_grid.month = ...
    repmat(permute(SOCATv2020_grid.month,[3 2 1]),...
    size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),1);
SOCATv2020_grid.latitude = repmat(SOCATv2020_grid.lat,1,1,max(size(SOCATv2020_grid.month_since_1998)));
SOCATv2020_grid.longitude = repmat(SOCATv2020_grid.lon,1,1,max(size(SOCATv2020_grid.month_since_1998)));

%% Determine distance to shore
disp('Determining distance from shore');
SOCATv2020_grid.distance_from_shore = ...
    dist2coast(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1));
SOCATv2020_grid.distance_from_shore = ...
    repmat(SOCATv2020_grid.distance_from_shore,1,1,max(size(SOCATv2020_grid.month_since_1998)));

%% Create indices for land, coastal product, and open-ocean product
SOCATv2020_grid.mask_land  = landmask(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360);

%% Obtain sea surface salinity from ECCO reanalysis
% Import ECCO2 SSS
disp('Obtaining ECCO2 sea surface salinity');
path = '/Volumes/2TB Hard Drive/ECCO2_SSS/';
load(strcat(path,'ECCO2_SSS_monthly_averaged_data.mat')); clear path
% Format latitude and longitude
ECCO_SSS.latitude = repmat(ECCO_SSS.lat',size(ECCO_SSS.sss_mon,1),1,size(ECCO_SSS.sss_mon,3));
ECCO_SSS.longitude = repmat(ECCO_SSS.lon,1,size(ECCO_SSS.sss_mon,2),size(ECCO_SSS.sss_mon,3));
% Cut down dataset to CCS limits
lonidx = ECCO_SSS.lon >= lonmin & ECCO_SSS.lon <= lonmax;
latidx = ECCO_SSS.lat >= latmin & ECCO_SSS.lat <= latmax;
ECCO_SSS.sss_mon = ECCO_SSS.sss_mon(lonidx,latidx,:);
ECCO_SSS.latitude = ECCO_SSS.latitude(lonidx,latidx,:);
ECCO_SSS.longitude = ECCO_SSS.longitude(lonidx,latidx,:);
% Match time frame of SOCAT data
ECCO_SSS.month_since_1998 = (ECCO_SSS.year_mon-1998).*12 + ECCO_SSS.mon_mon;
idx = ECCO_SSS.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      ECCO_SSS.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
ECCO_SSS.sss_mon = ECCO_SSS.sss_mon(:,:,idx);
ECCO_SSS.latitude = ECCO_SSS.latitude(:,:,idx);
ECCO_SSS.longitude = ECCO_SSS.longitude(:,:,idx);
% Save into SOCAT gridded structure (no lat-lon re-configuration necessary)
SOCATv2020_grid.SSS = double(ECCO_SSS.sss_mon);
% Interpolate over some gaps in SSS dataset
for t = 1:max(SOCATv2020_grid.month_since_1998)
    idx = ~isnan(SOCATv2020_grid.SSS(:,:,t)) & ~SOCATv2020_grid.mask_land;
    lon_tmp = SOCATv2020_grid.longitude(:,:,t);
    lat_tmp = SOCATv2020_grid.latitude(:,:,t);
    sss_tmp = SOCATv2020_grid.SSS(:,:,t);
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),sss_tmp(idx));
    idx = isnan(SOCATv2020_grid.SSS(:,:,t)) & ~SOCATv2020_grid.mask_land;
    sss_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    SOCATv2020_grid.SSS(:,:,t) = sss_tmp;
end
% Eliminate lake values
SOCATv2020_grid.SSS(SOCATv2020_grid.latitude > 40 & ...
    SOCATv2020_grid.longitude(:,:,1) > 240) = NaN;

% figure; pcolor(ECCO_SSS.longitude(:,:,1),ECCO_SSS.latitude(:,:,1),ECCO_SSS.sss_mon(:,:,1)); colorbar; caxis([30 36]);
% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.SSS(:,:,1)); colorbar; caxis([30 36]);

clear lonidx latidx idx interp lat_tmp lon_tmp sss_tmp t

%% Obtain sea surface height from ECCO reanalysis
% Import ECCO2 SSH
disp('Obtaining ECCO2 sea surface height');
path = '/Volumes/2TB Hard Drive/ECCO2_SSH/';
load(strcat(path,'ECCO2_SSH_monthly_averaged_data.mat')); clear path
% Format latitude and longitude
ECCO_SSH.latitude = repmat(ECCO_SSH.lat',size(ECCO_SSH.ssh_mon,1),1,size(ECCO_SSH.ssh_mon,3));
ECCO_SSH.longitude = repmat(ECCO_SSH.lon,1,size(ECCO_SSH.ssh_mon,2),size(ECCO_SSH.ssh_mon,3));
% Cut down dataset to CCS limits
lonidx = ECCO_SSH.lon >= lonmin & ECCO_SSH.lon <= lonmax;
latidx = ECCO_SSH.lat >= latmin & ECCO_SSH.lat <= latmax;
ECCO_SSH.ssh_mon = ECCO_SSH.ssh_mon(lonidx,latidx,:);
ECCO_SSH.latitude = ECCO_SSH.latitude(lonidx,latidx,:);
ECCO_SSH.longitude = ECCO_SSH.longitude(lonidx,latidx,:);
% Match time frame of SOCAT data
ECCO_SSH.month_since_1998 = (ECCO_SSH.year_mon-1998).*12 + ECCO_SSH.mon_mon;
idx = ECCO_SSH.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      ECCO_SSH.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
ECCO_SSH.ssh_mon = ECCO_SSH.ssh_mon(:,:,idx);
ECCO_SSH.latitude = ECCO_SSH.latitude(:,:,idx);
ECCO_SSH.longitude = ECCO_SSH.longitude(:,:,idx);
% Save into SOCAT gridded structure (no lat-lon re-configuration necessary)
SOCATv2020_grid.SSH = double(ECCO_SSH.ssh_mon);
% Interpolate over some gaps in SSH dataset (2-D, lat and lon)
for t = 1:max(SOCATv2020_grid.month_since_1998)
    idx = ~isnan(SOCATv2020_grid.SSH(:,:,t)) & ~SOCATv2020_grid.mask_land;
    lon_tmp = SOCATv2020_grid.longitude(:,:,t);
    lat_tmp = SOCATv2020_grid.latitude(:,:,t);
    ssh_tmp = SOCATv2020_grid.SSH(:,:,t);
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),ssh_tmp(idx));
    idx = isnan(SOCATv2020_grid.SSH(:,:,t)) & ~SOCATv2020_grid.mask_land;
    ssh_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    SOCATv2020_grid.SSH(:,:,t) = ssh_tmp;
end
% Eliminate lake values
SOCATv2020_grid.SSH(SOCATv2020_grid.latitude > 40 & ...
    SOCATv2020_grid.longitude(:,:,1) > 240) = NaN;

% figure; pcolor(ECCO_SSH.longitude(:,:,1),ECCO_SSH.latitude(:,:,1),ECCO_SSH.ssh_mon(:,:,1)); colorbar; caxis([30 36]);
% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.SSH(:,:,1)); colorbar; caxis([30 36]);

clear lonidx latidx idx interp lat_tmp lon_tmp ssh_tmp t

%% Obtain sea surface temperature from OISSTv2
% Import OISSTv2
disp('Obtaining OISSTv2 sea surface temperature');
path = '/Volumes/2TB Hard Drive/OISSTv2/';
load(strcat(path,'OISSTv2_monthly_averaged_data.mat')); clear path
% Format latitude and longitude
OISST.latitude = repmat(OISST.lat',size(OISST.sst_mon,1),1,size(OISST.sst_mon,3));
OISST.longitude = repmat(OISST.lon,1,size(OISST.sst_mon,2),size(OISST.sst_mon,3));
% Cut down dataset to CCS limits
lonidx = OISST.lon >= lonmin & OISST.lon <= lonmax;
latidx = OISST.lat >= latmin & OISST.lat <= latmax;
OISST.sst_mon = OISST.sst_mon(lonidx,latidx,:);
OISST.latitude = OISST.latitude(lonidx,latidx,:);
OISST.longitude = OISST.longitude(lonidx,latidx,:);
% Match time frame of SOCAT data
OISST.month_since_1998 = (OISST.year_mon-1998).*12 + OISST.mon_mon;
idx = OISST.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      OISST.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
OISST.sst_mon = OISST.sst_mon(:,:,idx);
OISST.latitude = OISST.latitude(:,:,idx);
OISST.longitude = OISST.longitude(:,:,idx);
% Save into SOCAT gridded structure (no lat-lon re-configuration necessary)
SOCATv2020_grid.SST = OISST.sst_mon;
% Interpolate over some gaps in SST dataset (2-D, lat and lon)
for t = 1:max(SOCATv2020_grid.month_since_1998)
    idx = ~isnan(SOCATv2020_grid.SST(:,:,t));% & ~SOCATv2020_grid.mask_land;
    lon_tmp = SOCATv2020_grid.longitude(:,:,t);
    lat_tmp = SOCATv2020_grid.latitude(:,:,t);
    sst_tmp = SOCATv2020_grid.SST(:,:,t);
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),sst_tmp(idx));
    idx = isnan(SOCATv2020_grid.SST(:,:,t));% & ~SOCATv2020_grid.mask_land;
    sst_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    SOCATv2020_grid.SST(:,:,t) = sst_tmp;
end
% Eliminate lake values
SOCATv2020_grid.SST(SOCATv2020_grid.latitude > 40 & ...
    SOCATv2020_grid.longitude(:,:,1) > 240) = NaN;

%figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.SST(:,:,1)); colorbar; caxis([5 30]);

clear path latidx lonidx t idx interp

%% Obtain sea surface chlorophyll from satellite measurements
disp('Obtaining SeaWiFS and MODIS surface chlorophyll');
path = '/Volumes/2TB Hard Drive/SATELLITE_DATA/';
load(strcat(path,'CHL.mat')); clear path
CHL.latitude = repmat(CHL.lat',1,size(CHL.chl,2),size(CHL.chl,3));
CHL.longitude = repmat(CHL.lon,size(CHL.chl,1),1,size(CHL.chl,3));
% Cut down dataset to CCS limits
lonidx = CHL.lon >= lonmin-360 & CHL.lon <= lonmax-360;
latidx = CHL.lat >= latmin & CHL.lat <= latmax;
CHL.chl = CHL.chl(latidx,lonidx,:);
CHL.latitude = CHL.latitude(latidx,lonidx,:);
CHL.longitude = CHL.longitude(latidx,lonidx,:);
CHL.longitude = CHL.longitude + 360; % Convert to 360 degrees longitude
% Match time frame of SOCAT data
CHL.date = datevec(CHL.time);
CHL.month_since_1998 = (CHL.date(:,1)-1998).*12 + CHL.date(:,2);
idx = CHL.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      CHL.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
CHL.chl = CHL.chl(:,:,idx);
CHL.latitude = CHL.latitude(:,:,idx);
CHL.longitude = CHL.longitude(:,:,idx);
disp('Interpolating surface chlorophyll to SOCAT grid');
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2020_grid.month_since_1998)
    interp = griddedInterpolant(flipud(CHL.longitude(:,:,t))',flipud(CHL.latitude(:,:,t))',flipud(CHL.chl(:,:,t))');
    SOCATv2020_grid.CHL(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
end
% Eliminate lake values
SOCATv2020_grid.CHL(SOCATv2020_grid.latitude > 40 & ...
    SOCATv2020_grid.longitude(:,:,1) > 240) = NaN;
% Interpolate over some gaps in CHL dataset (1-D, time)
for g = 1:size(SOCATv2020_grid.CHL,1)
    for h = 1:size(SOCATv2020_grid.CHL,2)
        if sum(~isnan(SOCATv2020_grid.CHL(g,h,:))) > 1
        Chl = squeeze(SOCATv2020_grid.CHL(g,h,:));
        idx = ~isnan(Chl);
        Chlfit = interp1(SOCATv2020_grid.month_since_1998(idx),...
            Chl(idx),SOCATv2020_grid.month_since_1998,'pchip','extrap');
        SOCATv2020_grid.CHL(g,h,:) = Chlfit;
        end
    end
end
% Replace values less than zero with value close to zero
SOCATv2020_grid.CHL(SOCATv2020_grid.CHL<0) = 0.001;

% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.CHL(:,:,1)); colorbar; caxis([0 2]);

clear latidx lonidx t idx interp Chl Chlfit

%% Obtain wind speed from ERA5 re-analysis
disp('Obtaining ERA5 re-analysis winds');
importERA5
ERA5.latitude = repmat(ERA5.lat',size(ERA5.speed,1),1,size(ERA5.speed,3));
ERA5.longitude = repmat(ERA5.lon+360,1,size(ERA5.speed,2),size(ERA5.speed,3));
% Match time frame of SOCAT data
ERA5.date = datevec(ERA5.date);
ERA5.month_since_1998 = (ERA5.date(:,1)-1998).*12 + ERA5.date(:,2);
idx = ERA5.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      ERA5.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
ERA5.speed = ERA5.speed(:,:,idx);
ERA5.u10 = ERA5.u10(:,:,idx);
ERA5.v10 = ERA5.v10(:,:,idx);
ERA5.latitude = ERA5.latitude(:,:,idx);
ERA5.longitude = ERA5.longitude(:,:,idx);
disp('Interpolating wind speed to SOCAT grid');
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2020_grid.month_since_1998)
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.speed(:,:,t)));
    SOCATv2020_grid.wind_speed(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.u10(:,:,t)));
    SOCATv2020_grid.u10(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
    interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.v10(:,:,t)));
    SOCATv2020_grid.v10(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
end

% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.wind_speed(:,:,1)); colorbar; caxis([0 12]);

clear t idx interp

% %% Obtain wind speed from CCMP
% disp('Obtaining CCMP re-analysis winds');
% importCCMP
% ERA5.latitude = repmat(ERA5.lat',size(ERA5.speed,1),1,size(ERA5.speed,3));
% ERA5.longitude = repmat(ERA5.lon+360,1,size(ERA5.speed,2),size(ERA5.speed,3));
% % Match time frame of SOCAT data
% ERA5.date = datevec(ERA5.date);
% ERA5.month_since_1998 = (ERA5.date(:,1)-1998).*12 + ERA5.date(:,2);
% idx = ERA5.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
%       ERA5.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
% ERA5.speed = ERA5.speed(:,:,idx);
% ERA5.u10 = ERA5.u10(:,:,idx);
% ERA5.v10 = ERA5.v10(:,:,idx);
% ERA5.latitude = ERA5.latitude(:,:,idx);
% ERA5.longitude = ERA5.longitude(:,:,idx);
% disp('Interpolating wind speed to SOCAT grid');
% % Interpolate onto SOCAT grid
% for t = 1:max(SOCATv2020_grid.month_since_1998)
%     interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.speed(:,:,t)));
%     SOCATv2020_grid.wind_speed(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
%     interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.u10(:,:,t)));
%     SOCATv2020_grid.u10(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
%     interp = griddedInterpolant(fliplr(ERA5.longitude(:,:,t)),fliplr(ERA5.latitude(:,:,t)),fliplr(ERA5.v10(:,:,t)));
%     SOCATv2020_grid.v10(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
% end
% 
% % figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.wind_speed(:,:,1)); colorbar; caxis([0 12]);
% 
% clear t idx interp

%% Obtain bathymetry from ETOPO2
disp('Obtaining ETOPO2 bathymetry');
importETOPO2
ETOPO2.latitude = repmat(ETOPO2.lat',size(ETOPO2.bottomdepth,1),1);
ETOPO2.longitude = repmat(ETOPO2.lon,1,size(ETOPO2.bottomdepth,2));
% Cut down dataset to CCS limits
lonidx = ETOPO2.lon >= lonmin-360 & ETOPO2.lon <= lonmax-360;
latidx = ETOPO2.lat >= latmin & ETOPO2.lat <= latmax;
ETOPO2.bottomdepth = ETOPO2.bottomdepth(lonidx,latidx,:);
ETOPO2.latitude = ETOPO2.latitude(lonidx,latidx,:);
ETOPO2.longitude = ETOPO2.longitude(lonidx,latidx,:);
ETOPO2.longitude = ETOPO2.longitude + 360; % Convert to 360 degrees longitude
% Interpolate onto SOCAT grid
disp('Interpolating bathymetry to SOCAT grid');
interp = griddedInterpolant(ETOPO2.longitude,ETOPO2.latitude,ETOPO2.bottomdepth);
SOCATv2020_grid.bottomdepth = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
% Extend to time frame of SOCAT data
SOCATv2020_grid.bottomdepth = repmat(SOCATv2020_grid.bottomdepth,1,1,max(size(SOCATv2020_grid.month_since_1998)));

% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.bottomdepth(:,:,1)); colorbar;

clear path latidx lonidx interp

%% Obtain mixed layer depth from HYCOM model
disp('Obtaining HYCOM mixed layer depth');
importMLD
MLD.latitude = repmat(MLD.lat',1,size(MLD.mld,2),size(MLD.mld,3));
MLD.longitude = repmat(MLD.lon,size(MLD.mld,1),1,size(MLD.mld,3));
% Cut down dataset to CCS limits
lonidx = MLD.lon >= lonmin-360 & MLD.lon <= lonmax-360;
latidx = MLD.lat >= latmin & MLD.lat <= latmax;
MLD.mld = MLD.mld(latidx,lonidx,:);
MLD.latitude = MLD.latitude(latidx,lonidx,:);
MLD.longitude = MLD.longitude(latidx,lonidx,:);
MLD.longitude = MLD.longitude + 360; % Convert to 360 degrees longitude
% Match time frame of SOCAT data
MLD.date = datevec(MLD.time);
MLD.month_since_1998 = (MLD.date(:,1)-1998).*12 + MLD.date(:,2);
idx = MLD.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
      MLD.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
MLD.mld = MLD.mld(:,:,idx);
MLD.latitude = MLD.latitude(:,:,idx);
MLD.longitude = MLD.longitude(:,:,idx);
% Interpolate onto SOCAT grid
disp('Interpolating mixed layer depth to SOCAT grid');
for t = 1:max(SOCATv2020_grid.month_since_1998)
    interp = griddedInterpolant(flipud(MLD.longitude(:,:,t))',flipud(MLD.latitude(:,:,t))',flipud(MLD.mld(:,:,t))');
    SOCATv2020_grid.MLD(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
end
% Interpolate over some gaps in MLD dataset (2-D, lat and lon)
for t = 1:max(SOCATv2020_grid.month_since_1998)
    idx = ~isnan(SOCATv2020_grid.MLD(:,:,t)) & ~SOCATv2020_grid.mask_land;
    lon_tmp = SOCATv2020_grid.longitude(:,:,t);
    lat_tmp = SOCATv2020_grid.latitude(:,:,t);
    mld_tmp = SOCATv2020_grid.MLD(:,:,t);
    interp = scatteredInterpolant(lon_tmp(idx),lat_tmp(idx),mld_tmp(idx));
    idx = isnan(SOCATv2020_grid.MLD(:,:,t)) & ~SOCATv2020_grid.mask_land;
    mld_tmp(idx) = interp(lon_tmp(idx),lat_tmp(idx));
    SOCATv2020_grid.MLD(:,:,t) = mld_tmp;
end
% Eliminate lake values
SOCATv2020_grid.MLD(SOCATv2020_grid.latitude > 40 & ...
    SOCATv2020_grid.longitude(:,:,1) > 240) = NaN;
% Replace values less than zero with value close to zero
SOCATv2020_grid.MLD(SOCATv2020_grid.MLD<0) = 0.001;

% figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.MLD(:,:,1)); colorbar; caxis([0 150]);

clear latidx lonidx t idx interp

% %% Obtain atmospheric pCO2 from Jena CarboScope
% disp('Obtaining Jena CarboScope atmospheric pCO2');
% importCarboScope
% CarboScope.latitude = repmat(CarboScope.lat',size(CarboScope.paCO2,1),1,size(CarboScope.paCO2,3));
% CarboScope.longitude = repmat(CarboScope.lon,1,size(CarboScope.paCO2,2),size(CarboScope.paCO2,3));
% % Cut down dataset to CCS limits
% lonidx = CarboScope.lon >= lonmin-360 & CarboScope.lon <= lonmax-360;
% latidx = CarboScope.lat >= latmin & CarboScope.lat <= latmax;
% CarboScope.paCO2 = CarboScope.paCO2(lonidx,latidx,:);
% CarboScope.latitude = CarboScope.latitude(lonidx,latidx,:);
% CarboScope.longitude = CarboScope.longitude(lonidx,latidx,:);
% CarboScope.longitude = CarboScope.longitude + 360; % Convert to 360 degrees longitude
% CarboScope.latitude = CarboScope.latitude(:,:,1);
% CarboScope.longitude = CarboScope.longitude(:,:,1);
% % Match time frame of SOCAT data
% CarboScope.date = datevec(CarboScope.date);
% CarboScope.month_since_1998 = (CarboScope.date(:,1)-1998).*12 + CarboScope.date(:,2);
% idx = CarboScope.month_since_1998 >= min(SOCATv2020_grid.month_since_1998) & ...
%       CarboScope.month_since_1998 <= max(SOCATv2020_grid.month_since_1998);
% CarboScope.date = CarboScope.date(idx,1:3);
% CarboScope.month_since_1998 = CarboScope.month_since_1998(idx);
% CarboScope.paCO2 = CarboScope.paCO2(:,:,idx);
% % Determine monthly averages
% CarboScope.paCO2_mon = ...
%     nan(size(CarboScope.longitude,1),size(CarboScope.latitude,2),max(SOCATv2020_grid.month_since_1998));
% for t = 1:max(SOCATv2020_grid.month_since_1998)
%     idx = CarboScope.month_since_1998 == t;
%     CarboScope.paCO2_mon(:,:,t) = mean(CarboScope.paCO2(:,:,idx),3,'omitnan');
% end
% % Interpolate onto SOCAT grid
% disp('Interpolating atmospheric pCO2 to SOCAT grid');
% for t = 1:max(SOCATv2020_grid.month_since_1998)
%     interp = griddedInterpolant(CarboScope.longitude(:,:),CarboScope.latitude(:,:),CarboScope.paCO2_mon(:,:,t));
%     SOCATv2020_grid.paCO2(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
% end
% SOCATv2020_grid.paCO2 = double(SOCATv2020_grid.paCO2);
% 
% % figure; pcolor(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.paCO2(:,:,1)); colorbar; %caxis([0 150]);
% 
% clear latidx lonidx t idx interp

%% Obtain atmospheric pressure from NCEP
disp('Obtaining NCEP atmospheric pressure');
path = '/Volumes/2TB Hard Drive/SATELLITE_DATA/NCEP_PRES/'; % this can be changed as necessary
% Read netcdf file
NCEP = netcdfreader(strcat(path,'mslp.mon.mean.nc'));
NCEP.latitude = repmat(NCEP.lat',size(NCEP.mslp,1),1);
NCEP.longitude = repmat(NCEP.lon,1,size(NCEP.mslp,2));
% Calculate date
NCEP.date = datevec(datenum([repmat(1800,size(NCEP.time,1),1) ones(size(NCEP.time,1),1) ones(size(NCEP.time,1),1) ...
    NCEP.time zeros(size(NCEP.time,1),1) zeros(size(NCEP.time,1),1)]));
NCEP.months_since_1998 = (NCEP.date(:,1) - 1998).*12 + NCEP.date(:,2);
% Cut down dataset to CCS limits
lonidx = NCEP.lon >= lonmin & NCEP.lon <= lonmax;
latidx = NCEP.lat >= latmin & NCEP.lat <= latmax;
NCEP.mslp = NCEP.mslp(lonidx,latidx,:);
NCEP.latitude = NCEP.latitude(lonidx,latidx);
NCEP.longitude = NCEP.longitude(lonidx,latidx);
% Index to 1998-2019
index = NCEP.months_since_1998 > 0 & NCEP.months_since_1998 <= max(SOCATv2020_grid.month_since_1998);
NCEP.months_since_1998 = NCEP.months_since_1998(index);
NCEP.mslp = NCEP.mslp(:,:,index);
% Interpolate onto SOCAT grid
disp('Interpolating atmospheric to SOCAT grid');
for t = 1:max(SOCATv2020_grid.month_since_1998)
    interp = griddedInterpolant(fliplr(NCEP.longitude(:,:)),fliplr(NCEP.latitude(:,:)),fliplr(NCEP.mslp(:,:,t)));
    SOCATv2020_grid.mslp(:,:,t) = interp(SOCATv2020_grid.longitude(:,:,1),SOCATv2020_grid.latitude(:,:,1));
end
% Convert pascals to atmospheres
SOCATv2020_grid.mslp = double(SOCATv2020_grid.mslp)./101325;

clear t interp lonidx latidx

%% Obtain atmospheric pCO2 from NOAA MBL product
disp('Obtaining NOAA MBL atmospheric pCO2');
% Open and scan file
file = fopen('/Users/Sharp/Documents/DATA/MBL_1998_2019.txt');
NOAA_MBL = textscan(file,'%f','Delimiter',',','CommentStyle','#');
file = fclose(file);
NOAA_MBL = reshape(cell2mat(NOAA_MBL),83,[]);
% Process data within file
MBL.year = repmat(NOAA_MBL(1,:),(size(NOAA_MBL,1)-1)./2,1);
MBL.CO2 = NOAA_MBL(2:2:end,:);
MBL.err = NOAA_MBL(3:2:end,:);
% Define latitudes
latsin = [-1.00  -0.95  -0.90  -0.85  -0.80  -0.75  -0.70  -0.65  -0.60 ...
          -0.55  -0.50  -0.45  -0.40  -0.35  -0.30  -0.25  -0.20  -0.15 ...
          -0.10  -0.05   0.00   0.05   0.10   0.15   0.20   0.25   0.30 ...
           0.35   0.40   0.45   0.50   0.55   0.60   0.65   0.70   0.75 ...
           0.80   0.85   0.90   0.95   1.00];
MBL.lat = asind(latsin)';
% Define year fraction interval based on monthly slices
MBL.year_mon = repmat(min(min(MBL.year)):1/12:max(max(MBL.year)),size(MBL.CO2,1),1);
MBL.year_mon = MBL.year_mon(:,1:end-1);
MBL.year_mon = MBL.year_mon + 1/24;
% Interpolate to monthly time slices 
MBL.CO2_mon = nan(size(MBL.year_mon));
MBL.err_mon = nan(size(MBL.year_mon));
for l = 1:size(MBL.CO2,1)
    MBL.CO2_mon(l,:) = interp1(MBL.year(l,:),MBL.CO2(l,:),MBL.year_mon(l,:));
    MBL.err_mon(l,:) = interp1(MBL.year(l,:),MBL.err(l,:),MBL.year_mon(l,:));
end
% Interpolate monthly values to gridded latitude interval
SOCATv2020_grid.xCO2_atm = nan(size(SOCATv2020_grid.lat,2),size(MBL.CO2_mon,2));
SOCATv2020_grid.xCO2_err_atm = nan(size(SOCATv2020_grid.lat,2),size(MBL.CO2_mon,2));
for m = 1:size(MBL.CO2_mon,2)
    SOCATv2020_grid.xCO2_atm(:,m) = interp1(MBL.lat,MBL.CO2_mon(:,m),SOCATv2020_grid.lat(1,:)');
    SOCATv2020_grid.xCO2_err_atm(:,m) = interp1(MBL.lat,MBL.err_mon(:,m),SOCATv2020_grid.lat(1,:)');
end
% Replicate across longitudes
SOCATv2020_grid.xCO2_atm = repmat(permute(SOCATv2020_grid.xCO2_atm,[3 1 2]),size(SOCATv2020_grid.lon,1),1,1);
SOCATv2020_grid.xCO2_err_atm = repmat(permute(SOCATv2020_grid.xCO2_err_atm,[3 1 2]),size(SOCATv2020_grid.lon,1),1,1);
% Calculate pCO2 from xCO2 (with vapor pressure correction)
SOCATv2020_grid.vapor_pressure = vpress(SOCATv2020_grid.SSS,SOCATv2020_grid.SST);
SOCATv2020_grid.pCO2_atm = SOCATv2020_grid.xCO2_atm.*...
    (SOCATv2020_grid.mslp-SOCATv2020_grid.vapor_pressure);

clear file l t z MBL* latsin 

%% Clean up

clear CarboScope CHL ECCO_SSS ECCO_SSH ERA5 ETOPO2 MLD OISST
