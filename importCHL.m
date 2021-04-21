% Import_CHL_SeaWiFS_MODIS
path1 = '/Volumes/2TB Hard Drive/SATELLITE_DATA/SeaWiFS_CHL/'; % this can be changed as necessary
path2 = '/Volumes/2TB Hard Drive/SATELLITE_DATA/MODIS_CHL/CHLs_MONTHLY_6TH/'; % this can be changed as necessary

% Define month and year limits
year = 1998:1:2020;
monthdaynorm = [1 32 60 91 121 152 182 213 244 274 305 335];
monthdayleap = [1 32 61 92 122 153 183 214 245 275 306 336];

% Import netcdf files for 1998-2002 (SeaWiFS)
for y=find(year==1998):find(year==2002)
    if year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
       year(y) == 2012 || year(y) == 2016 || year(y) == 2020
        monthday = monthdayleap;
    else
        monthday = monthdaynorm;
    end
    for m=1:max(size(monthday))
        try
        CHL.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',m))) = ...
            hdfreader(strcat(path1,'chl.',num2str(year(y)),sprintf('%03d',monthday(m)),'.hdf'));
        catch
        end
    end
end
% Clean up
clear m y path1

% Import netcdf files for 2003-2020 (MODIS)
for y=find(year==2003):find(year==2020)
    if year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
       year(y) == 2012 || year(y) == 2016 || year(y) == 2020
        monthday = monthdayleap;
    else
        monthday = monthdaynorm;
    end
    for m=1:max(size(monthday))
        try
        CHL.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',m))) = ...
            hdfreader(strcat(path2,'chl.',num2str(year(y)),sprintf('%03d',monthday(m)),'.hdf'));
        catch
        end
    end
end
% Clean up
clear m y path2 monthdayleap monthdaynorm

% Concatenate into 3-D varaibles
CHL.chl = nan(1080,2160,max(size(year))*12); % preallocate 
for y=1:max(size(year))
    for m=1:max(size(monthday))
        try
        CHL.chl(:,:,(y-1)*12+m) = ...
            CHL.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',m))).chl;
        catch
        end
    end
end
CHL.chl(CHL.chl==-9999)=NaN;
% Add latitude and longitude to structure
[CHL.lat,CHL.lon] = OSU_Lat_Lon(size(CHL.chl,1),size(CHL.chl,2));
% Add time variable to structure
year = repmat(year,12,1);
CHL.time = datenum([year(:) repmat([1:12]',size(year,2),1) ones(size(year,1)*size(year,2),1)]);

% Remove extraneous fields
fields = fieldnames(CHL);
idx = startsWith(fields,'y');
CHL = rmfield(CHL,fields(idx));
% Clean up
clear m y monthday year idx fields

