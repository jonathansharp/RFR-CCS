disp('Gridding SOCATv2020 observations');

%% Construct lat-lon grid
% Establish latitude and longitude minima and maxima
latmin = 15; latmax = 60;
lonmin = 220; lonmax = 255;
% Create grid
[lat,lon] = meshgrid(latmin+0.125:0.25:latmax,lonmin+0.125:0.25:lonmax);

dataset = fieldnames(SOCATv2020);
%dataset = {'all' 'train' 'test'};

for n = 1:numel(dataset)

%% Determine bin number of each data point
[cnt_lon,Xbins,Xnum] = histcounts(SOCATv2020.(dataset{n}).longitude,lonmin:0.25:lonmax);
[cnt_lat,Ybins,Ynum] = histcounts(SOCATv2020.(dataset{n}).latitude,latmin:0.25:latmax);
[cnt_mon,Zbins,Znum] = histcounts(SOCATv2020.(dataset{n}).month_since_1998,0.5:1:264.5);

%% Accumulate 3D grid by applying function to SOCAT values with bin numbers that match grid cells
subs = [Xnum, Ynum, Znum];
sz = [size(lon,1),size(lat,2),max(size(1:1:264))];
SOCATv2020_grid.(dataset{n}).count_ncruise = accumarray(subs, SOCATv2020.(dataset{n}).cruise, sz, @(x) numel(unique(x)), NaN);
SOCATv2020_grid.(dataset{n}).pco2_count_nobs = accumarray(subs, SOCATv2020.(dataset{n}).pCO2, sz, @numel, NaN);
SOCATv2020_grid.(dataset{n}).pco2_ave_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).pCO2, sz, @mean, NaN);
SOCATv2020_grid.(dataset{n}).pco2_std_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).pCO2, sz, @std, NaN);
SOCATv2020_grid.(dataset{n}).pco2_max_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).pCO2, sz, @max, NaN);
SOCATv2020_grid.(dataset{n}).pco2_min_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).pCO2, sz, @min, NaN);
SOCATv2020_grid.(dataset{n}).sst_count_nobs = accumarray(subs, SOCATv2020.(dataset{n}).temperature, sz, @numel, NaN);
SOCATv2020_grid.(dataset{n}).sst_ave_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).temperature, sz, @mean, NaN);
SOCATv2020_grid.(dataset{n}).sst_std_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).temperature, sz, @std, NaN);
SOCATv2020_grid.(dataset{n}).sst_max_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).temperature, sz, @max, NaN);
SOCATv2020_grid.(dataset{n}).sst_min_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).temperature, sz, @min, NaN);
SOCATv2020_grid.(dataset{n}).salinity_count_nobs = accumarray(subs, SOCATv2020.(dataset{n}).salinity, sz, @numel, NaN);
SOCATv2020_grid.(dataset{n}).salinity_ave_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).salinity, sz, @mean, NaN);
SOCATv2020_grid.(dataset{n}).salinity_std_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).salinity, sz, @std, NaN);
SOCATv2020_grid.(dataset{n}).salinity_max_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).salinity, sz, @max, NaN);
SOCATv2020_grid.(dataset{n}).salinity_min_unwtd = accumarray(subs, SOCATv2020.(dataset{n}).salinity, sz, @min, NaN);

%% Determine cruise-weighted means and standard deviations
% If more than one cruise is represented in a given grid cell, replace
% the unweighted value with a cruise-weighted value

% Pre-allocate with unweighted values
SOCATv2020_grid.(dataset{n}).pco2_ave_weighted  = SOCATv2020_grid.(dataset{n}).pco2_ave_unwtd;
SOCATv2020_grid.(dataset{n}).pco2_std_weighted  = SOCATv2020_grid.(dataset{n}).pco2_std_unwtd;
SOCATv2020_grid.(dataset{n}).sst_ave_weighted  = SOCATv2020_grid.(dataset{n}).sst_ave_unwtd;
SOCATv2020_grid.(dataset{n}).sst_std_weighted  = SOCATv2020_grid.(dataset{n}).sst_std_unwtd;
SOCATv2020_grid.(dataset{n}).salinity_ave_weighted  = SOCATv2020_grid.(dataset{n}).salinity_ave_unwtd;
SOCATv2020_grid.(dataset{n}).salinity_std_weighted  = SOCATv2020_grid.(dataset{n}).salinity_std_unwtd;
% Determine monthly means for individual years
for a = 1:size(lon,1)
    for b = 1:size(lat,2)
        for c = 1:264
            if SOCATv2020_grid.(dataset{n}).count_ncruise(a,b,c) > 1
                % index to specific grid cell
                idx = SOCATv2020.(dataset{n}).longitude >= lon(a,b) - 0.125 & ...
                      SOCATv2020.(dataset{n}).longitude < lon(a,b) + 0.125 & ...
                      SOCATv2020.(dataset{n}).latitude  >= lat(a,b) - 0.125 & ...
                      SOCATv2020.(dataset{n}).latitude  < lat(a,b) + 0.125 & ...
                      SOCATv2020.(dataset{n}).month_since_1998 == c;
                cruises = unique(SOCATv2020.(dataset{n}).expocode(idx));
                cruiselist = SOCATv2020.(dataset{n}).expocode(idx);
                pco2 = SOCATv2020.(dataset{n}).pCO2(idx);
                temperature = SOCATv2020.(dataset{n}).temperature(idx);
                salinity = SOCATv2020.(dataset{n}).salinity(idx);
                pco22 = nan(numel(cruises),1);
                temperature2 = nan(numel(cruises),1);
                salinity2 = nan(numel(cruises),1);
                for k=1:numel(cruises)
                    cruiseidx = strcmp(cruiselist,cruises(k));
                    pco22(k) = mean(pco2(cruiseidx));
                    temperature2(k) = mean(temperature(cruiseidx));
                    salinity2(k) = mean(salinity(cruiseidx));
                end
                SOCATv2020_grid.(dataset{n}).pco2_ave_weighted(a,b,c)  = mean(pco22);
                SOCATv2020_grid.(dataset{n}).pco2_std_weighted(a,b,c)  = std(pco22);
                SOCATv2020_grid.(dataset{n}).sst_ave_weighted(a,b,c)  = mean(temperature2);
                SOCATv2020_grid.(dataset{n}).sst_std_weighted(a,b,c)  = std(temperature2);
                SOCATv2020_grid.(dataset{n}).salinity_ave_weighted(a,b,c)  = mean(salinity2);
                SOCATv2020_grid.(dataset{n}).salinity_std_weighted(a,b,c)  = std(salinity2);
            end
        end
    end
end

end

SOCATv2020_grid.lat = lat; SOCATv2020_grid.lon = lon;
SOCATv2020_grid.month_since_1998 = 1:264;
SOCATv2020_grid.month_since_1998 = SOCATv2020_grid.month_since_1998';

clear a b c cnt* cruise* dataset idx k lat lon n pco* sal* subs sz temp* X* Y* Z*
