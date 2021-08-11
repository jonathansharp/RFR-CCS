%% Remove trend from time series using monthly domain mean pCO2
area_weights = SOCATv2021_grid.area_km2;
area_weights = repmat(area_weights,1,1,size(SOCATv2021_grid.all.pco2_ave_weighted,3));
area_weights(isnan(SOCATv2021_grid.all.pco2_ave_weighted)) = NaN;
% Calculate area-weighted domain mean:
SOCATv2021_grid.all.pco2_dom_mean = ...
    squeeze(sum(sum(SOCATv2021_grid.all.pco2_ave_weighted.*area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
% Calculate area-weighted domain standard deviation:
SOCATv2021_grid.all.pco2_dom_std = ...
    squeeze(sum(sum(SOCATv2021_grid.all.pco2_std_weighted.*area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
% Fit trend:
[yf,yr,x] = leastsq(SOCATv2021_grid.month_since_1998,SOCATv2021_grid.all.pco2_dom_mean,...
    0,0,0);
% Remove difference from mean for each month:
for m = 1:max(SOCATv2021_grid.month_since_1998)
SOCATv2021_grid.all.pco2_ave_weighted_detrend(:,:,m) = SOCATv2021_grid.all.pco2_ave_weighted(:,:,m) + (mean(yf,'omitnan') - yf(m,:));
end
% Calculate area-weighted detrended domain mean:
SOCATv2021_grid.all.pco2_dom_mean_detrend = ...
    squeeze(sum(sum(SOCATv2021_grid.all.pco2_ave_weighted_detrend.*area_weights,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));

%% Visualize domain mean over time
% figure; hold on;
% errorbar(SOCATv2021_grid.month_since_1998,...
%     SOCATv2021_grid.all.pco2_dom_mean,SOCATv2021_grid.all.pco2_dom_std,...
%     'LineStyle','none');
% scatter(SOCATv2021_grid.month_since_1998,...
%     SOCATv2021_grid.all.pco2_dom_mean,20,'k','filled');
% plot(SOCATv2021_grid.month_since_1998,yf,'k','LineWidth',2);
% xlabel('Month since Jan. 1 1998');
% ylabel('Mean Domain Surface pCO2'); hold off;

%% Visualize de-trended domain mean over time
% figure; hold on;
% errorbar(SOCATv2021_grid.month_since_1998,...
%     SOCATv2021_grid.all.pco2_dom_mean_detrend,SOCATv2021_grid.all.pco2_dom_std,...
%     'LineStyle','none');
% scatter(SOCATv2021_grid.month_since_1998,...
%     SOCATv2021_grid.all.pco2_dom_mean_detrend,20,'k','filled');
% plot(SOCATv2021_grid.month_since_1998,repmat(mean(yf,'omitnan'),max(size(yf)),1),'k','LineWidth',2);
% xlabel('Month since Jan. 1 1998');
% ylabel('Mean Domain Surface pCO2'); hold off;

%% Determine monthly means for one annual cycle
% Pre-allocate
SOCATv2021_grid.all.pco2_ave_weighted_clim  = ...
    nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
SOCATv2021_grid.all.pco2_std_weighted_clim   = ...
    nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
% Determine monthly means for one annual cycle
for m = 1:12
    SOCATv2021_grid.all.pco2_ave_weighted_clim(:,:,m) = ...
        mean(SOCATv2021_grid.all.pco2_ave_weighted_detrend(:,:,m:12:end),3,'omitnan');
    SOCATv2021_grid.all.pco2_std_weighted_clim(:,:,m) = ...
        mean(SOCATv2021_grid.all.pco2_std_weighted(:,:,m:12:end),3,'omitnan');
end

%% Determine monthly means for one annual cycle (up to 2015)
% Pre-allocate
SOCATv2021_grid.all.pco2_ave_weighted_clim_2015  = ...
    nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
SOCATv2021_grid.all.pco2_std_weighted_clim_2015   = ...
    nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
% Determine monthly means for one annual cycle
for m = 1:12
    SOCATv2021_grid.all.pco2_ave_weighted_clim_2015(:,:,m) = ...
        mean(SOCATv2021_grid.all.pco2_ave_weighted_detrend(:,:,m:12:(2015-1998+1).*12),3,'omitnan');
    SOCATv2021_grid.all.pco2_std_weighted_clim_2015(:,:,m) = ...
        mean(SOCATv2021_grid.all.pco2_std_weighted(:,:,m:12:(2015-1998+1).*12),3,'omitnan');
end

%% Clean up workspace
clear m x yf yr area_weights
