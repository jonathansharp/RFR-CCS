%% Create and save NetCDF of RFR-CCS

fname = ['RFR-CCS_v1_' strrep(strrep(datestr(clock),':','-'),' ','-')];

% Create file and variables
nccreate([fname '.nc'],'pco2','Format','classic','Dimensions', ...
    {'longitude',size(SOCATv2021_grid.pco2_RF,1),...
     'latitude',size(SOCATv2021_grid.pco2_RF,2),...
      'time',size(SOCATv2021_grid.pco2_RF,3)},...
      'Datatype','single');
nccreate([fname '.nc'],'longitude','Dimensions',{'longitude'},'Datatype','single');
nccreate([fname '.nc'],'latitude','Dimensions',{'latitude'},'Datatype','single');
nccreate([fname '.nc'],'time','Dimensions',{'time'},'Datatype','single');
nccreate([fname '.nc'],'sst','Dimensions',{'longitude','latitude','time'},'Datatype','single');
nccreate([fname '.nc'],'sss','Dimensions',{'longitude','latitude','time'},'Datatype','single');
% nccreate([fname '.nc'],'Fco2','Dimensions',{'longitude','latitude','time'},'Datatype','single');

% Write variables
ncwrite([fname '.nc'],'pco2',SOCATv2021_grid.pco2_RF);
ncwrite([fname '.nc'],'longitude',SOCATv2021_grid.lon(:,1,1)-360);
ncwrite([fname '.nc'],'latitude',SOCATv2021_grid.lat(1,:,1)');
ncwrite([fname '.nc'],'time',SOCATv2021_grid.month_since_1998-0.5);
ncwrite([fname '.nc'],'sst',SOCATv2021_grid.SST);
ncwrite([fname '.nc'],'sss',SOCATv2021_grid.SSS);
% ncwrite([fname '.nc'],'Fco2',SOCATv2021_grid.Fco2_RF_ERA5);

% Write attributes
ncwriteatt([fname '.nc'],'pco2','long_name','CO2 partial pressure at SST');
ncwriteatt([fname '.nc'],'pco2','units','uatm');
ncwriteatt([fname '.nc'],'pco2','Source','Sharp et al (in prep)');
ncwriteatt([fname '.nc'],'pco2','Missing_Value',NaN);
ncwriteatt([fname '.nc'],'longitude','units','degrees east');
ncwriteatt([fname '.nc'],'latitude','units','degrees north');
ncwriteatt([fname '.nc'],'time','units','months since 1998-1-1');
ncwriteatt([fname '.nc'],'sst','long_name','sea surface temperature');
ncwriteatt([fname '.nc'],'sst','units','degC');
ncwriteatt([fname '.nc'],'sst','Missing_Value',NaN);
ncwriteatt([fname '.nc'],'sst','Source','OISSTv2.1');
ncwriteatt([fname '.nc'],'sss','long_name','sea surface salinity');
ncwriteatt([fname '.nc'],'sss','units','N/A');
ncwriteatt([fname '.nc'],'sss','Missing_Value',NaN);
ncwriteatt([fname '.nc'],'sss','Source','ECCO2');
% ncwriteatt([fname '.nc'],'Fco2','long_name','Sea to air CO2 flux');
% ncwriteatt([fname '.nc'],'Fco2','units','mol C m^-2 month^-1');
% ncwriteatt([fname '.nc'],'Fco2','Missing_Value',NaN);
% ncwriteatt([fname '.nc'],'Fco2','Source','Sharp et al (in prep)');
% ncwriteatt([fname '.nc'],'Fco2','kw','Ho et al. (2006)');
% ncwriteatt([fname '.nc'],'Fco2','K0','Weiss (1974)');
% ncwriteatt([fname '.nc'],'Fco2','wind','ERA5');

%ncdisp([fname '.nc'])


%% Create and save Matlab file of RFR-CCS

pco2 = single(SOCATv2021_grid.pco2_RF);
pco2_units = 'uatm';
longitude = single(SOCATv2021_grid.lon(:,1,1)-360);
latitude = single(SOCATv2021_grid.lat(1,:,1)');
time = single(SOCATv2021_grid.month_since_1998-0.5);
time_units = 'months_since_1998-1-1';
sst = single(SOCATv2021_grid.SST);
sst_units = 'degC';
sss = single(SOCATv2021_grid.SSS);
sss_units = 'N/A';
Fco2 = SOCATv2021_grid.Fco2_RF_ERA5;
Fco2_units = 'mol_C_per_m2_per_month';
save([fname '.mat'],'pco2','pco2_units','longitude','latitude','time',...
    'time_units','sst','sst_units','sss','sss_units');

%matfile('RFR-CCS_v1.mat')