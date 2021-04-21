% ImportEN4
path = '/Volumes/2TB Hard Drive/EN4/'; % this can be changed as necessary

% Define month and year limits
year = 1998:1:2020;
month = 1:1:12;

% Import netcdf files
for y=1:max(size(year))
    for m=1:max(size(month))
        try
        EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))) = ...
            netcdfreader(strcat(path,'EN.4.2.1.f.analysis.l09.',num2str(year(y)),sprintf('%02d',month(m)),'.nc'));
        catch
        end
    end
end
% Clean up
clear m y path

% Select only surface layer
for y=1:max(size(year))
    for m=1:max(size(month))
        try
            vars = fieldnames(EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))));
            for v=1:max(size(vars)) 
                EN4temp = EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).(vars{v});
                if ndims(EN4temp) == 3
                    EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).(vars{v}) = EN4temp(:,:,1);
                end
            end
        catch
        end
    end
end
% Clean up
clear m y EN4temp v vars

% Concatenate into 3-D varaibles
EN4.SST        = nan(360,173,23*12-1); % preallocate 
EN4.SST_uncer  = nan(360,173,23*12-1); % preallocate 
EN4.SST_weight = nan(360,173,23*12-1); % preallocate 
EN4.SSS        = nan(360,173,23*12-1); % preallocate 
EN4.SSS_uncer  = nan(360,173,23*12-1); % preallocate 
EN4.SSS_weight = nan(360,173,23*12-1); % preallocate 
EN4.time       = nan(23*12-1,1); % preallocate 
for y=1:max(size(year))
    for m=1:max(size(month))
        try
        EN4.SST(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).temperature;
        EN4.SST_uncer(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).temperature_uncertainty;
        EN4.SST_weight(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).temperature_observation_weights;
        EN4.SSS(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).salinity;
        EN4.SSS_uncer(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).salinity_uncertainty;
        EN4.SSS_weight(:,:,(y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).salinity_observation_weights;
        EN4.time((y-1)*12+m) = ...
            EN4.(strcat('y',num2str(year(y)))).(strcat('m',sprintf('%02d',month(m)))).time;
        catch
        end
    end
end
% Add latitude and longitude to structure
EN4.lat = EN4.y2020.m01.lat;
EN4.lon = EN4.y2020.m01.lon;

% Remove extraneous fields
fields = fieldnames(EN4);
idx = startsWith(fields,'y');
EN4 = rmfield(EN4,fields(idx));

% Clean up
clear m y month year fields idx