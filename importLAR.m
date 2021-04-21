% This script converts the data contained in the netcdf file associated 
% to Laruelle et al. 2017 into a 216 by 720 by 1440 matrix.
% The script also converts the scaled data into uatm using a scaling factor

ncid = netcdf.open('/Users/sharp/Documents/DATA/Laurelle CO2/Coastal_SOM_FFN_2017.nc','NOWRITE');

Time=netcdf.getVar(ncid,0);
Lat=netcdf.getVar(ncid,1);
Lon=netcdf.getVar(ncid,2);
Scaling_factor=netcdf.getVar(ncid,3);  
Var=netcdf.getVar(ncid,4);  
netcdf.close(ncid)

Var=single(Var);    % converts data from integer into single precision real
Var=Var/Scaling_factor; % converts data back into uatm using scaling factor
N=Var==-999;	% identifies cells covered in ice and replaces data with NaN
Var(N)=NaN;

pCO2=NaN*ones(216,720,1440);
l=length(Lon);
t=length(Time);
lon=[-179.875:.25:179.875];
lat=[-89.875:.25:89.875];

for i=1:1:t         % creates 216 by 720 by 1440 matrix
    for j=1:1:l
        pCO2(i,Lat(j)*4+360.5,Lon(j)*4+720.5)=Var(j,i);
    end
end

% Create date variable
LAR.date = datenum([repmat(1998,t,1) [1:t]' ones(t,1)]);

% save Coastal_pCO2.mat pCO2 -v7.3

% figure
% pcolor(lon,lat,squeeze(nanmean(pCO2)));shading flat
% title('mean pCO2 over 1998-2015 in uatm')
% xlabel('Latitude')
% ylabel('Longitude')
% colorbar;caxis([230 450]);

% Save pCO2
LAR.pCO2 = pCO2;

% Average by month
LAR.pCO2_mon=nan(12,720,1440);
LAR.pCO2_mon_std=nan(12,720,1440);
for t = 1:12
    LAR.pCO2_mon(t,:,:)=mean(LAR.pCO2(t:12:end,:,:));
    LAR.pCO2_mon_std(t,:,:)=std(LAR.pCO2(t:12:end,:,:));
end

% Adjust dimensions
LAR.pCO2 = permute(LAR.pCO2,[3 2 1]);
LAR.pCO2_mon = permute(LAR.pCO2_mon,[3 2 1]);
LAR.pCO2_mon_std = permute(LAR.pCO2_mon_std,[3 2 1]);

% Save lat-lon and wrap longitudes to 360 degrees
LAR.latitude = repmat(lat,size(lon,2),1);
LAR.longitude = repmat(lon',1,size(lat,2));
LAR.longitude(LAR.longitude<0) = LAR.longitude(LAR.longitude<0) + 360;

% Reorder to match Landschutzer longitude order
LAR.pCO2 = cat(1,LAR.pCO2(size(LAR.pCO2,1)/2+1:end,:,:),...
    LAR.pCO2(1:size(LAR.pCO2,1)/2,:,:));
LAR.pCO2_mon = cat(1,LAR.pCO2_mon(size(LAR.pCO2_mon,1)/2+1:end,:,:),...
    LAR.pCO2_mon(1:size(LAR.pCO2_mon,1)/2,:,:));
LAR.pCO2_mon_std = cat(1,LAR.pCO2_mon_std(size(LAR.pCO2_mon_std,1)/2+1:end,:,:),...
    LAR.pCO2_mon_std(1:size(LAR.pCO2_mon_std,1)/2,:,:));
LAR.longitude = cat(1,LAR.longitude(size(LAR.longitude,1)/2+1:end,:),...
    LAR.longitude(1:size(LAR.longitude,1)/2,:));

% Cut to CCS
% Latitude and longitude limits
latmin = 15; latmax = 60;
lonmin = 220; lonmax = 255;
% latitude and longitude indices
latidx = LAR.latitude(:,:) >= latmin & LAR.latitude(:,:) <= latmax;
lonidx = LAR.longitude(:,:) >= lonmin & LAR.longitude(:,:) <= lonmax;
% pCO2
idx = latidx & lonidx; idx = repmat(idx,1,1,size(LAR.pCO2,3));
LAR.pCO2 = reshape(LAR.pCO2(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),size(LAR.pCO2,3));
% Monthly pCO2
idx = latidx & lonidx; idx = repmat(idx,1,1,size(LAR.pCO2_mon,3));
LAR.pCO2_mon = reshape(LAR.pCO2_mon(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),size(LAR.pCO2_mon,3));
LAR.pCO2_mon_std = reshape(LAR.pCO2_mon_std(idx),max(max(sum(idx,1))),...
    max(max(sum(idx,2))),size(LAR.pCO2_mon_std,3));
% Latitude and longitude
idx = latidx & lonidx;
LAR.latitude = reshape(LAR.latitude(idx),max(sum(idx,1)),max(sum(idx,2)));
LAR.longitude = reshape(LAR.longitude(idx),max(sum(idx,1)),max(sum(idx,2)));

% Annual mean pCO2 and seasonal amplitude
LAR.pCO2_annmean = nanmean(LAR.pCO2_mon,3);
LAR.pCO2_amp = ...
max(LAR.pCO2_mon,[],3,'omitnan') - ...
min(LAR.pCO2_mon,[],3,'omitnan');

% Distance from shore
LAR.dist = dist2coast(LAR.latitude,LAR.longitude);
idx_dist = LAR.dist <= 100;
LAR.pCO2_dom_annmean = nanmean(nanmean(LAR.pCO2_annmean(idx_dist),1),2);
LAR.pCO2_dom_amp = nanmean(nanmean(LAR.pCO2_amp(idx_dist),1),2);

% % Plot Laruelle climatology
% figure; worldmap([round(latmin) latmax],[round(lonmin) round(lonmax)]);
% title('Annual Mean Laruelle pCO2')
% set(gcf,'Position',[617, 599, 420, 420])
% setm(gca,'ffacecolor',[0.94 0.97 1.0]);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
% contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
%     nanmean(LAR.pCO2_mon(:,:,:),3),...
%     300:2:500,'LineStyle','none');
% c=colorbar; colormap(jet); caxis([300 500]);
% c.Label.String = 'pCO2 (uatm)';

% % Plot Laruelle climatology amplitude
% figure; worldmap([round(latmin) latmax],[round(lonmin) round(lonmax)]);
% title('Laruelle pCO2 Amplitude')
% set(gcf,'Position',[617, 599, 420, 420])
% setm(gca,'ffacecolor',[0.94 0.97 1.0]);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
% contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
%     max(LAR.pCO2_mon(:,:,:),[],3,'omitnan')-...
%     min(LAR.pCO2_mon(:,:,:),[],3,'omitnan'),...
%     0:1:100,'LineStyle','none');
% c=colorbar; colormap(jet); caxis([0 100]);
% c.Label.String = 'pCO2 (uatm)';

clear i j l lat Lat lon Lon N ncid pCO2 Scaling_factor t Time Var latidx lonidx idx idx_dist
