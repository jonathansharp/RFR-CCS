disp('Gridding SOCATv2021 observations');

%% Construct lat-lon grid
% Establish latitude and longitude minima and maxima
latmin = round(min(SOCATv2021.all.latitude),0);
latmax = round(max(SOCATv2021.all.latitude),0);
lonmin = round(min(SOCATv2021.all.longitude),0);
lonmax = round(max(SOCATv2021.all.longitude),0);
% Create grid
[lat,lon] = meshgrid(latmin+0.125:0.25:latmax,lonmin+0.125:0.25:lonmax);

dataset = fieldnames(SOCATv2021);
%dataset = {'all' 'train' 'test'};

for n = 1:numel(dataset)

%% Determine bin number of each data point
[cnt_lon,Xbins,Xnum] = histcounts(SOCATv2021.(dataset{n}).longitude,lonmin:0.25:lonmax);
[cnt_lat,Ybins,Ynum] = histcounts(SOCATv2021.(dataset{n}).latitude,latmin:0.25:latmax);
[cnt_mon,Zbins,Znum] = histcounts(SOCATv2021.(dataset{n}).month_since_1998,0.5:1:276.5);

%% Accumulate 3D grid by applying function to SOCAT values with bin numbers that match grid cells
subs = [Xnum, Ynum, Znum];
sz = [size(lon,1),size(lat,2),max(size(1:1:276))];
SOCATv2021_grid.(dataset{n}).count_ncruise = accumarray(subs, SOCATv2021.(dataset{n}).cruise, sz, @(x) numel(unique(x)), NaN);
SOCATv2021_grid.(dataset{n}).pco2_count_nobs = accumarray(subs, SOCATv2021.(dataset{n}).pCO2, sz, @numel, NaN);
SOCATv2021_grid.(dataset{n}).pco2_ave_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).pCO2, sz, @mean, NaN);
SOCATv2021_grid.(dataset{n}).pco2_std_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).pCO2, sz, @std, NaN);
SOCATv2021_grid.(dataset{n}).pco2_max_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).pCO2, sz, @max, NaN);
SOCATv2021_grid.(dataset{n}).pco2_min_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).pCO2, sz, @min, NaN);
SOCATv2021_grid.(dataset{n}).sst_count_nobs = accumarray(subs, SOCATv2021.(dataset{n}).temperature, sz, @numel, NaN);
SOCATv2021_grid.(dataset{n}).sst_ave_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).temperature, sz, @mean, NaN);
SOCATv2021_grid.(dataset{n}).sst_std_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).temperature, sz, @std, NaN);
SOCATv2021_grid.(dataset{n}).sst_max_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).temperature, sz, @max, NaN);
SOCATv2021_grid.(dataset{n}).sst_min_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).temperature, sz, @min, NaN);
SOCATv2021_grid.(dataset{n}).salinity_count_nobs = accumarray(subs, SOCATv2021.(dataset{n}).salinity, sz, @numel, NaN);
SOCATv2021_grid.(dataset{n}).salinity_ave_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).salinity, sz, @mean, NaN);
SOCATv2021_grid.(dataset{n}).salinity_std_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).salinity, sz, @std, NaN);
SOCATv2021_grid.(dataset{n}).salinity_max_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).salinity, sz, @max, NaN);
SOCATv2021_grid.(dataset{n}).salinity_min_unwtd = accumarray(subs, SOCATv2021.(dataset{n}).salinity, sz, @min, NaN);

%% Determine cruise-weighted means and standard deviations
% If more than one cruise is represented in a given grid cell, replace
% the unweighted value with a cruise-weighted value

% Pre-allocate with unweighted values
SOCATv2021_grid.(dataset{n}).pco2_ave_weighted  = SOCATv2021_grid.(dataset{n}).pco2_ave_unwtd;
SOCATv2021_grid.(dataset{n}).pco2_std_weighted  = SOCATv2021_grid.(dataset{n}).pco2_std_unwtd;
SOCATv2021_grid.(dataset{n}).sst_ave_weighted  = SOCATv2021_grid.(dataset{n}).sst_ave_unwtd;
SOCATv2021_grid.(dataset{n}).sst_std_weighted  = SOCATv2021_grid.(dataset{n}).sst_std_unwtd;
SOCATv2021_grid.(dataset{n}).salinity_ave_weighted  = SOCATv2021_grid.(dataset{n}).salinity_ave_unwtd;
SOCATv2021_grid.(dataset{n}).salinity_std_weighted  = SOCATv2021_grid.(dataset{n}).salinity_std_unwtd;
SOCATv2021_grid.(dataset{n}).pco2_grid_uncert  = nan(size(SOCATv2021_grid.(dataset{n}).pco2_std_unwtd));
% Determine monthly means for individual years
for a = 1:size(lon,1)
    for b = 1:size(lat,2)
        for c = 1:276
            if SOCATv2021_grid.(dataset{n}).count_ncruise(a,b,c) > 1
                % index to specific grid cell
                idx = SOCATv2021.(dataset{n}).longitude >= lon(a,b) - 0.125 & ...
                      SOCATv2021.(dataset{n}).longitude < lon(a,b) + 0.125 & ...
                      SOCATv2021.(dataset{n}).latitude  >= lat(a,b) - 0.125 & ...
                      SOCATv2021.(dataset{n}).latitude  < lat(a,b) + 0.125 & ...
                      SOCATv2021.(dataset{n}).month_since_1998 == c;
                cruises = unique(SOCATv2021.(dataset{n}).expocode(idx));
                cruiselist = SOCATv2021.(dataset{n}).expocode(idx);
                pco2 = SOCATv2021.(dataset{n}).pCO2(idx);
                temperature = SOCATv2021.(dataset{n}).temperature(idx);
                salinity = SOCATv2021.(dataset{n}).salinity(idx);
                pco22 = nan(numel(cruises),1);
                temperature2 = nan(numel(cruises),1);
                salinity2 = nan(numel(cruises),1);
                for k=1:numel(cruises)
                    cruiseidx = strcmp(cruiselist,cruises(k));
                    pco22(k) = mean(pco2(cruiseidx));
                    temperature2(k) = mean(temperature(cruiseidx));
                    salinity2(k) = mean(salinity(cruiseidx));
                end
                SOCATv2021_grid.(dataset{n}).pco2_ave_weighted(a,b,c)  = mean(pco22);
                SOCATv2021_grid.(dataset{n}).pco2_std_weighted(a,b,c)  = std(pco22);
                SOCATv2021_grid.(dataset{n}).pco2_grid_uncert(a,b,c)  = std(pco22);
                SOCATv2021_grid.(dataset{n}).sst_ave_weighted(a,b,c)  = mean(temperature2);
                SOCATv2021_grid.(dataset{n}).sst_std_weighted(a,b,c)  = std(temperature2);
                SOCATv2021_grid.(dataset{n}).salinity_ave_weighted(a,b,c)  = mean(salinity2);
                SOCATv2021_grid.(dataset{n}).salinity_std_weighted(a,b,c)  = std(salinity2);
            end
        end
    end
end

end

SOCATv2021_grid.lat = lat; SOCATv2021_grid.lon = lon;
SOCATv2021_grid.month_since_1998 = 1:276;
SOCATv2021_grid.month_since_1998 = SOCATv2021_grid.month_since_1998';

clear a b c cnt* cruise* idx k lat lon n pco* sal* subs sz temp* X* Y* Z*

%% Determine area of each grid cell
disp('Determining area of each grid cell');
SOCATv2021_grid.area_km2 = ...
(((SOCATv2021_grid.lat + 0.125) - ...
    (SOCATv2021_grid.lat - 0.125)) .* 110.574) .* ... % latitude distance
(((SOCATv2021_grid.lon + 0.125) - ...
    (SOCATv2021_grid.lon - 0.125)) .* ...
    111.320.*cosd(SOCATv2021_grid.lat)); % longitude distance
