% ImportERA5hourly
path = '/Volumes/2TB Hard Drive/SATELLITE_DATA/ERA5/'; % this can be changed as necessary

% Read netcdf file
ERA5u = netcdfreader(strcat(path,'ERA5_CCS_3hourly_uwind.nc'));
ERA5v = netcdfreader(strcat(path,'ERA5_CCS_3hourly_vwind.nc'));

% Extract components
ERA5h.u10 = ERA5u.u10(:,:,1,:);
ERA5h.u10 = squeeze(permute(ERA5h.u10,[1 2 4 3]));
ERA5h.v10 = ERA5v.v10(:,:,1,:);
ERA5h.v10 = squeeze(permute(ERA5h.v10,[1 2 4 3]));

% Calculate wind speed
ERA5h.speed = sqrt(ERA5h.u10.^2 + ERA5h.v10.^2);

% Convert lat and lon to double
ERA5h.lat = double(ERA5u.latitude);
ERA5h.lon = double(ERA5u.longitude);

% Add date to structure
ERA5h.date = datenum([repmat(1900,max(size(ERA5u.time)),1) ...
             ones(max(size(ERA5u.time)),1) ...
             ones(max(size(ERA5u.time)),1) ...
             double(ERA5u.time) ,...
             zeros(max(size(ERA5u.time)),1) ...
             zeros(max(size(ERA5u.time)),1)]);

% Clean up
clear path ERA5u ERA5v