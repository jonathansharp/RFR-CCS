%% Uncertainty in pCO2 and CO2 Flux

% Uncertainty = u(measurement)^2 + u(grid)^2 + u(model)^2

%% Measurement uncertainty
% Measurement uncertainty is either 2 uatm or 5 uatm

% Convert dataset flags to numeric:
SOCATv2020.all.flagn = SOCATv2020.all.flag;
SOCATv2020.all.flagn(strcmp(SOCATv2020.all.flag,'A')) = {1};
SOCATv2020.all.flagn(strcmp(SOCATv2020.all.flag,'B')) = {2};
SOCATv2020.all.flagn(strcmp(SOCATv2020.all.flag,'C')) = {3};
SOCATv2020.all.flagn(strcmp(SOCATv2020.all.flag,'D')) = {4};
SOCATv2020.all.flagn = cell2mat(SOCATv2020.all.flagn);

% Percent with estimated uncertainty of 2 uatm vs. 5 uatm:
SOCATv2020.all.per_2uatm = ...
    sum(SOCATv2020.all.flagn == 1 | SOCATv2020.all.flagn == 2)./ ...
    size(SOCATv2020.all.flagn,1);
SOCATv2020.all.per_5uatm = 1 - SOCATv2020.all.per_2uatm;

% Average uncertainty:
u_meas = 2*SOCATv2020.all.per_2uatm + 5*SOCATv2020.all.per_5uatm

%% Grid uncertainty
% This is attributable to one value in a grid cell for a month not being
% representative of any location within that grid cell during any time
% within that month

% Determine bin number of each data point from original data
[cnt_lon,Xbins,Xnum] = histcounts(SOCATv2020.all.longitude,lonmin:0.25:lonmax);
[cnt_lat,Ybins,Ynum] = histcounts(SOCATv2020.all.latitude,latmin:0.25:latmax);
[cnt_mon,Zbins,Znum] = histcounts(SOCATv2020.all.month_since_1998,0.5:1:264.5);

% Index to coast
coast_idx = SOCATv2020_grid.distance_from_shore < 400;
SOCATv2020_grid.all.pco2_count_nobs_coast = ...
    SOCATv2020_grid.all.pco2_count_nobs(coast_idx);
SOCATv2020_grid.all.pco2_count_nobs_open = ...
    SOCATv2020_grid.all.pco2_count_nobs(~coast_idx);

% Take most densely populated grid cells (coastal)
[~,dense_coast] = maxk(SOCATv2020_grid.all.pco2_count_nobs_coast(:),10);

% Calculate standard deviation between all values in month (coastal)
pco2_std_coast = SOCATv2020_grid.all.pco2_std_unwtd(coast_idx);
mean(pco2_std_coast(dense_coast))

% Take most densely populated grid cells (open)
[~,dense_open] = maxk(SOCATv2020_grid.all.pco2_count_nobs_open(:),10);

% Calculate standard deviation between all values in month (open)
pco2_std_open = SOCATv2020_grid.all.pco2_std_unwtd(~coast_idx);
mean(pco2_std_open(dense_open))

% Grid pCO2 observations to coarse grid
gridSOCATv2020_coarse;

% Replicate coarse grid at resolution of fine grid
d = 1:size(SOCATv2020_grid_coarse.all.pco2_ave_weighted,1);
d = repelem(d,2);
e = 1:size(SOCATv2020_grid_coarse.all.pco2_ave_weighted,2);
e = repelem(e,2);
f = 1:size(SOCATv2020_grid_coarse.all.pco2_ave_weighted,3);
f = repelem(f,2);
for a = 1:size(SOCATv2020_grid.all.pco2_ave_weighted,1)
    for b = 1:size(SOCATv2020_grid.all.pco2_ave_weighted,2)
        for c = 1:size(SOCATv2020_grid.all.pco2_ave_weighted,3)
            SOCATv2020_grid_coarse.expanded.pco2_ave_weighted(a,b,c) = ...
                SOCATv2020_grid_coarse.all.pco2_ave_weighted(d(a),e(b),f(c));
        end
    end
end

% Compare coarse grid values to fine grid values (coast)
diff_coast = ...
    SOCATv2020_grid_coarse.expanded.pco2_ave_weighted(coast_idx) - ...
    SOCATv2020_grid.all.pco2_ave_weighted(coast_idx);
% Average difference between 0.5 x 0.5 2-month grid and 0.25 x 0.25 1-month grid
mean(diff_coast,'omitnan')
% Scale to ...

% Compare coarse grid values to fine grid values (open)
diff_open = ...
    SOCATv2020_grid_coarse.expanded.pco2_ave_weighted(~coast_idx) - ...
    SOCATv2020_grid.all.pco2_ave_weighted(~coast_idx);
% Average difference between 0.5 x 0.5 2-month grid and 0.25 x 0.25 1-month grid
mean(diff_open,'omitnan')
% Scale to ...

%% Model uncertainty
% This is attributable to the mismatch metween values predicted by the
% random forest regression model and gridded values

% Compare predicted values to test values (coastal)


% Compare predicted values to test values (open)




