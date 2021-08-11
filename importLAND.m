%% Import Landschutzer Climatology
path = 'Data/MPI-ULB-SOM_FFN_clim.nc';
LAND.latitude = ncread(path,'lat');
LAND.latitude_bounds = ncread(path,'lat_bnds');
LAND.longitude = ncread(path,'lon');
LAND.longitude_bounds = ncread(path,'lon_bnds');
LAND.time = ncread(path,'time');
LAND.pCO2 = ncread(path,'pco2');
LAND.pCO2 = ncread(path,'pco2');
% Change -999 to NaN in Landschutzer climatology
LAND.pCO2 = double(LAND.pCO2);
LAND.pCO2(LAND.pCO2==-999) = NaN;
% Expand Landschutzer climatology variables to same size
LAND.latitude = double(repmat(LAND.latitude',size(LAND.pCO2,1),1,size(LAND.pCO2,3)));
LAND.longitude = double(repmat(LAND.longitude,1,size(LAND.pCO2,2),size(LAND.pCO2,3)));
clear path

%% Process Landschutzer latitude and longitude
% Wrap longitudes to 360 degrees
LAND.longitude(LAND.longitude<0) = LAND.longitude(LAND.longitude<0) + 360;
% Re-order Landschutzer variables based on 360 degrees longitude
for t = 1:size(LAND.longitude,3)
    for c = 1:size(LAND.longitude,2)
        [~,id]=sort(LAND.longitude(:,c,t),1);
        longitude = LAND.longitude(:,c,t); LAND.longitude(:,c,t) = longitude(id);
        latitude  = LAND.latitude(:,c,t);  LAND.latitude(:,c,t)  = latitude(id);
        pCO2      = LAND.pCO2(:,c,t);      LAND.pCO2(:,c,t)      = pCO2(id);
    end
end

% Cut to CCS
% Latitude and longitude limits
latmin = 15; latmax = 60;
lonmin = 220; lonmax = 255;
% latitude and longitude indices
latidx = LAND.latitude(:,:,1) >= latmin & LAND.latitude(:,:,1) <= latmax;
lonidx = LAND.longitude(:,:,1) >= lonmin & LAND.longitude(:,:,1) <= lonmax;
idx = latidx & lonidx; idx = repmat(idx,1,1,size(LAND.pCO2,3));
% Monthly pCO2
LAND.pCO2 = reshape(LAND.pCO2(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),max(max(sum(idx,3))));
% Latitude and longitude
LAND.latitude = reshape(LAND.latitude(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),max(max(sum(idx,3))));
LAND.longitude = reshape(LAND.longitude(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),max(max(sum(idx,3))));

% Annual mean pCO2 and seasonal amplitude
LAND.pCO2_annmean = nanmean(LAND.pCO2,3);
LAND.pCO2_amp = ...
max(LAND.pCO2,[],3,'omitnan') - ...
min(LAND.pCO2,[],3,'omitnan');

% Distance from shore
LAND.dist = dist2coast(LAND.latitude(:,:,1),LAND.longitude(:,:,1));
idx_dist = LAND.dist <= 100;
LAND.pCO2_coast_annmean = nanmean(nanmean(LAND.pCO2_annmean(idx_dist),1),2);
LAND.pCO2_coast_amp = nanmean(nanmean(LAND.pCO2_amp(idx_dist),1),2);

% % Plot Landschutzer climatology
% figure; worldmap([round(latmin) latmax],[round(lonmin) round(lonmax)]);
% title('Annual Mean Landschutzer pCO2')
% set(gcf,'Position',[617, 599, 420, 420])
% setm(gca,'ffacecolor',[0.94 0.97 1.0]);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
% s=geoshape(lat(idx_dist),lon(idx_dist),'polygon')
% geoshow(s)
% contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
%     LAND.pCO2_annmean,...
%     300:2:500,'LineStyle','none');
% c=colorbar; colormap(jet); caxis([300 500]);
% c.Label.String = 'pCO2 (uatm)';
% 
% % Plot Landschutzer climatology amplitude
% figure; worldmap([round(latmin) latmax],[round(lonmin) round(lonmax)]);
% title('Landschutzer pCO2 Amplitude')
% set(gcf,'Position',[617, 599, 420, 420])
% setm(gca,'ffacecolor',[0.94 0.97 1.0]);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
% contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
%     LAND.pCO2_amp,...
%     0:1:100,'LineStyle','none');
% c=colorbar; colormap(jet); caxis([0 100]);
% c.Label.String = 'pCO2 (uatm)';

% % Determine land mask
% LAND.longitude180 = LAND.longitude;
% LAND.longitude180(LAND.longitude180>180) = LAND.longitude180(LAND.longitude180>180)-360;
% mask = landmask(LAND.latitude(:,:,1),LAND.longitude180(:,:,1));

clear t c id latitude longitude pCO2 latidx lonidx idx
