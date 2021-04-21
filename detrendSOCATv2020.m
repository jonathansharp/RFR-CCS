%% Remove trend from time series using monthly domain mean pCO2
% Calculate domain mean:
SOCATv2020_grid.all.pco2_dom_mean = squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_ave_weighted,1),2));
% Fit trend:
[yf,yr,x] = leastsq(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.all.pco2_dom_mean,...
    0,0,0);
% Remove difference from mean for each month:
for m = 1:max(SOCATv2020_grid.month_since_1998)
SOCATv2020_grid.all.pco2_ave_weighted_detrend(:,:,m) = SOCATv2020_grid.all.pco2_ave_weighted(:,:,m) + (nanmean(yf) - yf(m,:));
end

%% Visualize domain mean over time
figure; hold on;
errorbar(SOCATv2020_grid.month_since_1998,...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_ave_weighted(:,:,:),1),2)),...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_std_weighted(:,:,:),1),2)),...
    'LineStyle','none');
scatter(SOCATv2020_grid.month_since_1998,...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_ave_weighted(:,:,:),1),2)),...
    20,'k','filled');
plot(SOCATv2020_grid.month_since_1998,yf,'k','LineWidth',2);
xlabel('Month since Jan. 1 1998');
ylabel('Mean Domain Surface pCO2'); hold off;

%% Visualize de-trended domain mean over time
figure; hold on;
errorbar(SOCATv2020_grid.month_since_1998,...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_ave_weighted_detrend(:,:,:),1),2)),...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_std_weighted(:,:,:),1),2)),...
    'LineStyle','none');
scatter(SOCATv2020_grid.month_since_1998,...
    squeeze(nanmean(nanmean(SOCATv2020_grid.all.pco2_ave_weighted_detrend(:,:,:),1),2)),...
    20,'k','filled');
plot(SOCATv2020_grid.month_since_1998,repmat(nanmean(yf),max(size(yf)),1),'k','LineWidth',2);
xlabel('Month since Jan. 1 1998');
ylabel('Mean Domain Surface pCO2'); hold off;

%% Determine monthly means for one annual cycle
% Pre-allocate
SOCATv2020_grid.all.pco2_ave_weighted_clim  = ...
    nan(size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),12);
SOCATv2020_grid.all.pco2_std_weighted_clim   = ...
    nan(size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),12);
% Determine monthly means for one annual cycle
for m = 1:12
    SOCATv2020_grid.all.pco2_ave_weighted_clim(:,:,m) = mean(SOCATv2020_grid.all.pco2_ave_weighted_detrend(:,:,m:12:end),3,'omitnan');
    SOCATv2020_grid.all.pco2_std_weighted_clim(:,:,m) = mean(SOCATv2020_grid.all.pco2_std_weighted(:,:,m:12:end),3,'omitnan');
end

%% Determine monthly means for one annual cycle (up to 2015)
% Pre-allocate
SOCATv2020_grid.all.pco2_ave_weighted_clim_2015  = ...
    nan(size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),12);
SOCATv2020_grid.all.pco2_std_weighted_clim_2015   = ...
    nan(size(SOCATv2020_grid.lon,1),size(SOCATv2020_grid.lat,2),12);
% Determine monthly means for one annual cycle
for m = 1:12
    SOCATv2020_grid.all.pco2_ave_weighted_clim_2015(:,:,m) = mean(SOCATv2020_grid.all.pco2_ave_weighted_detrend(:,:,m:12:(2015-1998+1).*12),3,'omitnan');
    SOCATv2020_grid.all.pco2_std_weighted_clim_2015(:,:,m) = mean(SOCATv2020_grid.all.pco2_std_weighted(:,:,m:12:(2015-1998+1).*12),3,'omitnan');
end

clear m x yf yr
