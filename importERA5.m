% ImportERA5
path = '/Volumes/2TB Hard Drive/SATELLITE_DATA/ERA5/'; % this can be changed as necessary

% Read netcdf file
ERA5 = netcdfreader(strcat(path,'ERA5_CCS.nc'));

% Extract wind speed
ERA5.speed = ERA5.si10(:,:,1,:);
ERA5.speed = squeeze(permute(ERA5.speed,[1 2 4 3]));
ERA5.u10 = ERA5.u10(:,:,1,:);
ERA5.u10 = squeeze(permute(ERA5.u10,[1 2 4 3]));
ERA5.v10 = ERA5.v10(:,:,1,:);
ERA5.v10 = squeeze(permute(ERA5.v10,[1 2 4 3]));

% Convert lat and lon to double
ERA5.lat = double(ERA5.latitude);
ERA5.lon = double(ERA5.longitude);

% Add date to structure
ERA5.date = datenum([repmat(1900,max(size(ERA5.time)),1) ...
             ones(max(size(ERA5.time)),1) ...
             ones(max(size(ERA5.time)),1) ...
             double(ERA5.time) ,...
             zeros(max(size(ERA5.time)),1) ...
             zeros(max(size(ERA5.time)),1)]);

% Clean up
clear path 
% Remove extraneous fields
ERA5 = rmfield(ERA5,{'longitude' 'latitude' 'expver' 'time' 'si10'});