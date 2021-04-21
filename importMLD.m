% importMLD
path = '/Volumes/2TB Hard Drive/MLD/MLDs_MONTHLY_6TH/'; % this can be changed as necessary

% Define month and year limits
year = 1998:1:2020;
monthdaynorm = [1 32 60 91 121 152 182 213 244 274 305 335];
monthdayleap = [1 32 61 92 122 153 183 214 245 275 306 336];

% Import netcdf files
for y=1:max(size(year))
    if year(y) == 2000 || year(y) == 2004 || year(y) == 2008 || ...
       year(y) == 2012 || year(y) == 2016 || year(y) == 2020
        monthday = monthdayleap;
    else
        monthday = monthdaynorm;
    end
    for m=1:max(size(monthday))
        %try
        MLD.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',m))) = ...
            hdfreader(strcat(path,'mld.',num2str(year(y)),sprintf('%03d',monthday(m)),'.hdf'));
        %catch
        %end
    end
end
% Clean up
clear m y path1 monthdayleap monthdaynorm

% Concatenate into 3-D varaibles
MLD.mld = nan(1080,2160,max(size(year))*12); % preallocate 
for y=1:max(size(year))
    for m=1:max(size(monthday))
        try
        MLD.mld(:,:,(y-1)*12+m) = ...
            MLD.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',m))).mld;
        catch
        end
    end
end
MLD.mld(MLD.mld==-9999)=NaN;
% Add latitude and longitude to structure
[MLD.lat,MLD.lon] = OSU_Lat_Lon(size(MLD.mld,1),size(MLD.mld,2));
% Add time variable to structure
year = repmat(year,12,1);
MLD.time = datenum([year(:) repmat([1:12]',size(year,2),1) ones(size(year,1)*size(year,2),1)]);

% Remove extraneous fields
fields = fieldnames(MLD);
idx = startsWith(fields,'y');
MLD = rmfield(MLD,fields(idx));
% Clean up
clear m y monthday year idx fields path



