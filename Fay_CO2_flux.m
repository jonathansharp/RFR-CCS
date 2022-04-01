% compare flux from RFR-CCS to Fay et al. (2021) product

%% import Fay et al. pCO2 products
Fay.filler = ncread('Data/SeaFlux_v2021.04_spco2_filler_1990-2019.nc','spco2_scaled_climatology');
Fay.CMEMS_FFNN = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','CMEMS_FFNN');
Fay.CSIR_ML6 = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','CSIR_ML6');
Fay.JENA_MLS = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','JENA_MLS');
Fay.JMA_MLR = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','JMA_MLR');
Fay.MPI_SOMFFN = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','MPI_SOMFFN');
Fay.NIES_FNN = ncread('Data/SeaFlux_v2021.04_spco2_SOCOM_unfilled_1982-2019.nc','NIES_FNN');

%% import Fay et al. FCO2 product
Fay.Fco2 = ncread('Data/SeaFlux_v2021.04_fgco2_all_winds_products.nc','fgco2');
Fay.lat = ncread('Data/SeaFlux_v2021.04_fgco2_all_winds_products.nc','lat');
Fay.lon = ncread('Data/SeaFlux_v2021.04_fgco2_all_winds_products.nc','lon');
Fay.time = ncread('Data/SeaFlux_v2021.04_fgco2_all_winds_products.nc','time');
Fay.area = ncread('SeaFlux_v2021.04_area_ocean.nc','area_ocean');

%% limit to region corresponding to RFR-CCS
if lonmin > 180; lonmin = lonmin-360; end
if lonmax > 180; lonmax = lonmax-360; end
lonidx = Fay.lon >= lonmin & Fay.lon <= lonmax;
latidx = Fay.lat >= latmin & Fay.lat <= latmax;
Fay.filler = Fay.filler(lonidx,latidx,:);
Fay.CMEMS_FFNN = Fay.CMEMS_FFNN(lonidx,latidx,:);
Fay.CSIR_ML6 = Fay.CSIR_ML6(lonidx,latidx,:);
Fay.JENA_MLS = Fay.JENA_MLS(lonidx,latidx,:);
Fay.JMA_MLR = Fay.JMA_MLR(lonidx,latidx,:);
Fay.MPI_SOMFFN = Fay.MPI_SOMFFN(lonidx,latidx,:);
Fay.NIES_FNN = Fay.NIES_FNN(lonidx,latidx,:);
Fay.Fco2 = Fay.Fco2(lonidx,latidx,:,:,:);
Fay.area = Fay.area(lonidx,latidx);
Fay.latitude = repmat(Fay.lat(latidx)',sum(lonidx),1);
Fay.longitude = repmat(Fay.lon(lonidx),1,sum(latidx));

%% percent of climatology-filled grid cells
Fay.CMEMS_FFNN = 100.*(squeeze(sum(sum(~isnan(Fay.CMEMS_FFNN),2),1))./mean(sum(sum(~isnan(Fay.filler),2),1)));
Fay.CSIR_ML6 = 100.*(squeeze(sum(sum(~isnan(Fay.CSIR_ML6),2),1))./mean(sum(sum(~isnan(Fay.filler(:,:,1)),2),1)));
Fay.JENA_MLS = 100.*(squeeze(sum(sum(~isnan(Fay.JENA_MLS),2),1))./mean(sum(sum(~isnan(Fay.filler(:,:,1)),2),1)));
Fay.JMA_MLR = 100.*(squeeze(sum(sum(~isnan(Fay.JMA_MLR),2),1))./mean(sum(sum(~isnan(Fay.filler(:,:,1)),2),1)));
Fay.MPI_SOMFFN = 100.*(squeeze(sum(sum(~isnan(Fay.MPI_SOMFFN),2),1))./mean(sum(sum(~isnan(Fay.filler(:,:,1)),2),1)));
Fay.NIES_FNN = 100.*(squeeze(sum(sum(~isnan(Fay.NIES_FNN),2),1))./mean(sum(sum(~isnan(Fay.filler(:,:,1)),2),1)));

Fay.CMEMS_FFNN = mean(Fay.CMEMS_FFNN(Fay.CMEMS_FFNN~=0));
Fay.CSIR_ML6 = mean(Fay.CSIR_ML6(Fay.CSIR_ML6~=0));
Fay.JENA_MLS = mean(Fay.JENA_MLS(Fay.JENA_MLS~=0));
Fay.JMA_MLR = mean(Fay.JMA_MLR(Fay.JMA_MLR~=0));
Fay.MPI_SOMFFN = mean(Fay.MPI_SOMFFN(Fay.MPI_SOMFFN~=0));
Fay.NIES_FNN = mean(Fay.NIES_FNN(Fay.NIES_FNN~=0));

Average_spatial_coverage = mean([Fay.CMEMS_FFNN,Fay.CSIR_ML6,Fay.JENA_MLS,...
                                 Fay.JMA_MLR,Fay.MPI_SOMFFN,Fay.NIES_FNN])

%% match time frame of RFR-CCS
% calculate month since Jan. 1 1998
Fay.date = datenum([repmat(1981,length(Fay.time),1) ...
    repmat(12,length(Fay.time),1) 15+double(Fay.time)]);
Fay.month_since_1998 = datevec(Fay.date);
Fay.month_since_1998 = ...
    (Fay.month_since_1998(:,1)-1998).*12 + Fay.month_since_1998(:,2);
idx = Fay.month_since_1998 > 0;
Fay.Fco2 = Fay.Fco2(:,:,idx,:,:);
Fay.month_since_1998 = Fay.month_since_1998(idx);

%% average FCO2 across pCO2 products and extract ERA5 wind results
Fay.Fco2_mean = mean(Fay.Fco2(:,:,:,2,:),5);
Fay.Fco2_std = std(Fay.Fco2(:,:,:,2,:),[],5);
Fay.Fco2_CMEMS_FFNN = Fay.Fco2(:,:,:,2,1);
Fay.Fco2_CSIR_ML6 = Fay.Fco2(:,:,:,2,2);
Fay.Fco2_JENA_MLS = Fay.Fco2(:,:,:,2,3);
Fay.Fco2_JMA_MLR = Fay.Fco2(:,:,:,2,4);
Fay.Fco2_MPI_SOMFFN = Fay.Fco2(:,:,:,2,5);
Fay.Fco2_NIES_FNN = Fay.Fco2(:,:,:,2,6);

%% determine coastal and open-ocean indices
Fay.dist = dist2coast(Fay.latitude,Fay.longitude)
Fay.open_idx = Fay.dist >= 100 & Fay.area > 0;
Fay.coast_idx = Fay.dist < 100 & Fay.area > 0;
Fay.both_idx = Fay.open_idx | Fay.coast_idx;

%% take regional spatially-averaged CO2 flux
[Fay.Fco2_dom_mean,Fay.Fco2_dom_cum] = domain_mean(Fay.Fco2_mean,Fay.area,Fay.both_idx);

[Fay.Fco2_dom_mean_open,Fay.Fco2_dom_cum_open] = domain_mean(Fay.Fco2_mean,Fay.area,Fay.open_idx);
[Fay.Fco2_dom_mean_plus_open,Fay.Fco2_dom_cum_plus_open] = domain_mean(Fay.Fco2_mean+Fay.Fco2_std,Fay.area,Fay.open_idx);
[Fay.Fco2_dom_mean_minus_open,Fay.Fco2_dom_cum_minus_open] = domain_mean(Fay.Fco2_mean-Fay.Fco2_std,Fay.area,Fay.open_idx);

[Fay.Fco2_dom_mean_coast,Fay.Fco2_dom_cum_coast] = domain_mean(Fay.Fco2_mean,Fay.area,Fay.coast_idx);
[Fay.Fco2_dom_mean_plus_coast,Fay.Fco2_dom_cum_plus_coast] = domain_mean(Fay.Fco2_mean+Fay.Fco2_std,Fay.area,Fay.coast_idx);
[Fay.Fco2_dom_mean_minus_coast,Fay.Fco2_dom_cum_minus_coast] = domain_mean(Fay.Fco2_mean-Fay.Fco2_std,Fay.area,Fay.coast_idx);

[Fay.Fco2_CMEMS_FFNN_dom_mean_open,Fay.Fco2_CMEMS_FFNN_dom_cum_open] = ...
    domain_mean(Fay.Fco2_CMEMS_FFNN,Fay.area,Fay.open_idx);
[Fay.Fco2_CSIR_ML6_dom_mean_open,Fay.Fco2_CSIR_ML6_dom_cum_open] = ...
    domain_mean(Fay.Fco2_CSIR_ML6,Fay.area,Fay.open_idx);
[Fay.Fco2_JENA_MLS_dom_mean_open,Fay.Fco2_JENA_MLS_dom_cum_open] = ...
    domain_mean(Fay.Fco2_JENA_MLS,Fay.area,Fay.open_idx);
[Fay.Fco2_JMA_MLR_dom_mean_open,Fay.Fco2_JMA_MLR_dom_cum_open] = ...
    domain_mean(Fay.Fco2_JMA_MLR,Fay.area,Fay.open_idx);
[Fay.Fco2_MPI_SOMFFN_dom_mean_open,Fay.Fco2_MPI_SOMFFN_dom_cum_open] = ...
    domain_mean(Fay.Fco2_MPI_SOMFFN,Fay.area,Fay.open_idx);
[Fay.Fco2_NIES_FNN_dom_mean_open,Fay.Fco2_NIES_FNN_dom_cum_open] = ...
    domain_mean(Fay.Fco2_NIES_FNN,Fay.area,Fay.open_idx);

[Fay.Fco2_CMEMS_FFNN_dom_mean_coast,Fay.Fco2_CMEMS_FFNN_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_CMEMS_FFNN,Fay.area,Fay.coast_idx);
[Fay.Fco2_CSIR_ML6_dom_mean_coast,Fay.Fco2_CSIR_ML6_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_CSIR_ML6,Fay.area,Fay.coast_idx);
[Fay.Fco2_JENA_MLS_dom_mean_coast,Fay.Fco2_JENA_MLS_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_JENA_MLS,Fay.area,Fay.coast_idx);
[Fay.Fco2_JMA_MLR_dom_mean_coast,Fay.Fco2_JMA_MLR_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_JMA_MLR,Fay.area,Fay.coast_idx);
[Fay.Fco2_MPI_SOMFFN_dom_mean_coast,Fay.Fco2_MPI_SOMFFN_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_MPI_SOMFFN,Fay.area,Fay.coast_idx);
[Fay.Fco2_NIES_FNN_dom_mean_coast,Fay.Fco2_NIES_FNN_dom_cum_coast] = ...
    domain_mean(Fay.Fco2_NIES_FNN,Fay.area,Fay.coast_idx);

%% calculate climatology
for m = 1:12
    
    Fay.Fco2_dom_cum_clim_open(m) = mean(Fay.Fco2_dom_cum_open(m:12:end));
    Fay.Fco2_dom_cum_plus_clim_open(m) = mean(Fay.Fco2_dom_cum_plus_open(m:12:end));
    Fay.Fco2_dom_cum_minus_clim_open(m) = mean(Fay.Fco2_dom_cum_minus_open(m:12:end));
    Fay.Fco2_dom_cum_clim_coast(m) = mean(Fay.Fco2_dom_cum_coast(m:12:end));
    Fay.Fco2_dom_cum_plus_clim_coast(m) = mean(Fay.Fco2_dom_cum_plus_coast(m:12:end));
    Fay.Fco2_dom_cum_minus_clim_coast(m) = mean(Fay.Fco2_dom_cum_minus_coast(m:12:end));
    
    Fay.Fco2_dom_mean_clim_open(m) = mean(Fay.Fco2_dom_mean_open(m:12:end));
    Fay.Fco2_dom_mean_plus_clim_open(m) = mean(Fay.Fco2_dom_mean_plus_open(m:12:end));
    Fay.Fco2_dom_mean_minus_clim_open(m) = mean(Fay.Fco2_dom_mean_minus_open(m:12:end));
    Fay.Fco2_dom_mean_clim_coast(m) = mean(Fay.Fco2_dom_mean_coast(m:12:end));
    Fay.Fco2_dom_mean_plus_clim_coast(m) = mean(Fay.Fco2_dom_mean_plus_coast(m:12:end));
    Fay.Fco2_dom_mean_minus_clim_coast(m) = mean(Fay.Fco2_dom_mean_minus_coast(m:12:end));
    
    Fay.Fco2_CMEMS_FFNN_dom_mean_clim_open(m) = mean(Fay.Fco2_CMEMS_FFNN_dom_mean_open(m:12:end));
    Fay.Fco2_CSIR_ML6_dom_mean_clim_open(m) = mean(Fay.Fco2_CSIR_ML6_dom_mean_open(m:12:end));
    Fay.Fco2_JENA_MLS_dom_mean_clim_open(m) = mean(Fay.Fco2_JENA_MLS_dom_mean_open(m:12:end));
    Fay.Fco2_JMA_MLR_dom_mean_clim_open(m) = mean(Fay.Fco2_JMA_MLR_dom_mean_open(m:12:end));
    Fay.Fco2_MPI_SOMFFN_dom_mean_clim_open(m) = mean(Fay.Fco2_MPI_SOMFFN_dom_mean_open(m:12:end));
    Fay.Fco2_NIES_FNN_dom_mean_clim_open(m) = mean(Fay.Fco2_NIES_FNN_dom_mean_open(m:12:end));

    Fay.Fco2_CMEMS_FFNN_dom_mean_clim_coast(m) = mean(Fay.Fco2_CMEMS_FFNN_dom_mean_coast(m:12:end));
    Fay.Fco2_CSIR_ML6_dom_mean_clim_coast(m) = mean(Fay.Fco2_CSIR_ML6_dom_mean_coast(m:12:end));
    Fay.Fco2_JENA_MLS_dom_mean_clim_coast(m) = mean(Fay.Fco2_JENA_MLS_dom_mean_coast(m:12:end));
    Fay.Fco2_JMA_MLR_dom_mean_clim_coast(m) = mean(Fay.Fco2_JMA_MLR_dom_mean_coast(m:12:end));
    Fay.Fco2_MPI_SOMFFN_dom_mean_clim_coast(m) = mean(Fay.Fco2_MPI_SOMFFN_dom_mean_coast(m:12:end));
    Fay.Fco2_NIES_FNN_dom_mean_clim_coast(m) = mean(Fay.Fco2_NIES_FNN_dom_mean_coast(m:12:end));
    
    SOCATv2021_grid.Fco2_RF_dom_cum_clim(m) = mean(SOCATv2021_grid.Fco2_RF_dom_cum(m:12:end));
    SOCATv2021_grid.Fco2_RF_dom_cum_clim_open(m) = mean(SOCATv2021_grid.Fco2_RF_dom_cum_open(m:12:end));
    SOCATv2021_grid.Fco2_RF_dom_cum_clim_coast(m) = mean(SOCATv2021_grid.Fco2_RF_dom_cum_coast(m:12:end));
    
    SOCATv2021_grid.Fco2_RF_dom_mean_clim(m) = mean(SOCATv2021_grid.Fco2_RF_dom_mean(m:12:end));
    SOCATv2021_grid.Fco2_RF_dom_mean_clim_open(m) = mean(SOCATv2021_grid.Fco2_RF_dom_mean_open(m:12:end));
    SOCATv2021_grid.Fco2_RF_dom_mean_clim_coast(m) = mean(SOCATv2021_grid.Fco2_RF_dom_mean_coast(m:12:end));

end

%% produce figure (open-ocean climatology)
% Cumulative
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 0.3 0.5]);
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_cum_plus_clim_open./12 ...
    fliplr(Fay.Fco2_dom_cum_minus_clim_open./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
plot(1:12,Fay.Fco2_dom_cum_clim_open./12,'linewidth',6,'color',[0 0 0]);
plot(1:12,SOCATv2021_grid.Fco2_RF_dom_cum_clim_open./12,'linewidth',6,'color',[0 0 0.9]);
xlim([0 13]); xticks([1:12]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('Cumulative {\itF}CO_{2} (Tg C yr^{-1})','fontsize',22);
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_RFR_CCS_vs_Fay.jpg');
% Per area
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 0.3 0.5]);
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_open./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_open./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
plot(1:12,Fay.Fco2_dom_mean_clim_open./12,'linewidth',6,'color',[0 0 0]);
plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_open./1000,'linewidth',6,'color',[0 0 0.9]);
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('{\itF}CO_{2} (mol C m^{-2} yr^{-1})','fontsize',22);
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_RFR_CCS_vs_Fay.jpg');

%% produce figure (coastal climatology)
% Cumulative
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 0.3 0.5]);
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_cum_plus_clim_coast./12 ...
    fliplr(Fay.Fco2_dom_cum_minus_clim_coast./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
plot(1:12,Fay.Fco2_dom_cum_clim_coast./12,'linewidth',6,'color',[0 0 0]);
plot(1:12,SOCATv2021_grid.Fco2_RF_dom_cum_clim_coast./12,'linewidth',6,'color',[0 0 0.9]);
xlim([0 13]); xticks([1:12]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('Cumulative {\itF}CO_{2} (Tg C yr^{-1})','fontsize',22);
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_RFR_CCS_vs_Fay.jpg');
% Per area
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 0.3 0.5]);
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_coast./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_coast./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
plot(1:12,Fay.Fco2_dom_mean_clim_coast./12,'linewidth',6,'color',[0 0 0]);
plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_coast./1000,'linewidth',6,'color',[0 0 0.9]);
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('{\itF}CO_{2} (mol C m^{-2} yr^{-1})','fontsize',22);
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_RFR_CCS_vs_Fay.jpg');

%% produce figure (time series)
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 1 0.5]);
set(gca,'fontsize',18);
plot(Fay.month_since_1998,Fay.Fco2_dom_cum./12,'linewidth',3);
plot(SOCATv2021_grid.month_since_1998(1:264),...
    SOCATv2021_grid.Fco2_RF_dom_cum(1:264)./12,'linewidth',3);
xlim([-5 271]); xticks(2*12:5*12:23*12);
xticklabels({'2000' '2005' '2010' '2015' '2020'});
ylabel('Cumulative {\itF}CO_{2} (Tg C yr^{-1})','fontsize',26);
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_RFR_CCS_vs_Fay.jpg');

%% produce figure (spatial flux)
% figure; hold on;
% set(gcf,'units','normalized','outerposition',[0.0 0.5 0.5 0.5]);
% worldmap(latlims,lonlims);
% setm(gca,'ffacecolor',ocncol);
% setm(gca,'fontsize',fntsz);
% set(gca,'fontsize',fntsz);
% land = shaperead('landareas', 'UseGeoCoords', true);
% contourfm(Fay.latitude,Fay.longitude,...
%     mean(Fay.Fco2_mean,3,'omitnan'),...
%     -5:0.5:5,'LineStyle','none');
% geoshow(land, 'FaceColor',lndcol);
% caxis([-5 5]);
% colormap(cmocean('balance',20,'pivot',0));
% 
%% produce figure (spatial flux)
% figure; hold on;
% set(gcf,'units','normalized','outerposition',[0.0 0.5 0.5 0.5]);
% worldmap(latlims,lonlims);
% setm(gca,'ffacecolor',ocncol);
% setm(gca,'fontsize',fntsz);
% set(gca,'fontsize',fntsz);
% land = shaperead('landareas', 'UseGeoCoords', true);
% contourfm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
%     mean(SOCATv2021_grid.Fco2_RF_ERA5./1000.*12,3,'omitnan'),...
%     -5:0.5:5,'LineStyle','none');
% geoshow(land, 'FaceColor',lndcol);
% caxis([-5 5]);
% colormap(cmocean('balance',20,'pivot',0));

%% calculate annual mean flux
Fay_annual_mean = mean(Fay.Fco2_dom_cum)
Fay_annual_mean_open = mean(Fay.Fco2_dom_cum_open)
Fay_annual_mean_coast = mean(Fay.Fco2_dom_cum_coast)
%sum(sum(Fay.area))
RFR_CCS_annual_mean = mean(SOCATv2021_grid.Fco2_RF_dom_cum)
RFR_CCS_annual_mean_open = mean(SOCATv2021_grid.Fco2_RF_dom_cum_open)
RFR_CCS_annual_mean_coast = mean(SOCATv2021_grid.Fco2_RF_dom_cum_coast)
%sum(sum(SOCATv2021_grid.area_km2.*SOCATv2021_grid.percent_sea.*10^6))

%% embedded function

function [fco2_dom_mean,fco2_dom_cum] = domain_mean(fco2,area,idx)

area(~idx) = 0;

fco2_dom_mean = ... % mol C m^-2 yr^-1
    squeeze(sum(sum(fco2.*area,1,'omitnan'),2,'omitnan'))./...
    squeeze(sum(sum(area,1,'omitnan'),2,'omitnan'));

fco2_dom_cum = ... % Tg C yr^-1
    squeeze(sum(sum(fco2.*area,1,'omitnan'),2,'omitnan')).*...
    12.011./... % g/mol
    (10^12); % g/Tg

end