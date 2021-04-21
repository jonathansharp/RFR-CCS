%% Figures
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
%ocncol = [0.94 0.97 1.0];
ocncol = [1 1 1];
fntsz = 18;
%lndcol = [0.2 0.2 0.2];
lndcol = [1 1 1];

%% Plot on grid

% Plot annual mean gridded pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.lat,SOCATv2020_grid.lon,...
    nanmean(SOCATv2020_grid.all.pco2_ave_weighted_clim,3),...
    300:1:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
colormap(cmocean('haline'));
caxis([300 440]);
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_grid.jpg');

% Plot gridded pCO2 amplitude
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.lat,SOCATv2020_grid.lon,...
    max(SOCATv2020_grid.all.pco2_ave_weighted_clim(:,:,:),[],3,'omitnan')-...
    min(SOCATv2020_grid.all.pco2_ave_weighted_clim(:,:,:),[],3,'omitnan'),...
    0:1:100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([0 100]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_grid_amplitude.jpg');

%  Plot percentage of months represented in grid
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.lat,SOCATv2020_grid.lon,...
    100.*(sum(~isnan(SOCATv2020_grid.all.pco2_ave_weighted),3)./size(SOCATv2020_grid.all.pco2_ave_weighted,3)),...
    0:0.1:10,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([0 10]);
c.Label.String = '% of months';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_grid_month_percent.jpg');

%% Plot merged Landschutzer estimates
importLAND

% Plot annual mean Landschutzer pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
    nanmean(LAND.pCO2(:,:,:),3),...
    280:1:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([280 440]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_LAND.jpg');

% Plot Landschutzer pCO2 amplitude
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
    max(LAND.pCO2,[],3,'omitnan')-min(LAND.pCO2,[],3,'omitnan'),...
    0:1:100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(cmocean('thermal')); caxis([0 100]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_LAND_amplitude.jpg');

% Plot Landschutzer values vs SOCAT grid cells
figure;
set(gca,'fontsize',fntsz);
Index = ~isnan(SOCATv2020_grid.all.pco2_ave_weighted_clim_2015) & ~isnan(LAND.pCO2);
scatter(SOCATv2020_grid.all.pco2_ave_weighted_clim_2015(Index),...
    LAND.pCO2(Index),10,log10(SOCATv2020_grid.distance_from_shore(Index)),'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCAT gridded monthly mean pCO_{2} (\muatm)','fontsize',16);
xlim([100 1000]);
ylabel('Landschutzer et al. (2020) climatological pCO_{2} (\muatm)','fontsize',16);
ylim([100 1000]);
c=colorbar; colormap(cmocean('balance'));
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/Land_vs_SOCAT_scatter.jpg');

% Plot Landscutzer-predicted CO2 flux
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '/Users/sharp/Desktop/CO2_Flux_LAND.gif';
MonthNames={'January' 'February' 'March' 'April' 'May' 'June' 'July' 'August' 'September' 'October' 'November' 'December'};
for n = 1:12

worldmap(latlims,lonlims);
title(MonthNames(n),'fontsize',18);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
    mean(LAND.Fco2_ERA5(:,:,n:12:end),3,'omitnan'),...
    -2:0.05:2,'LineStyle','none');
LAND.longitude_condensed = repmat([220:1:255]',1,46);
LAND.latitude_condensed = repmat(15:1:60,36,1);
for i=1:216
    LAND.v10_condensed(:,:,i) = interp2(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
        SOCATv2020_grid.v10(:,:,i),...
        LAND.latitude_condensed,LAND.longitude_condensed);
    LAND.u10_condensed(:,:,i) = interp2(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
        SOCATv2020_grid.u10(:,:,i),...
        LAND.latitude_condensed,LAND.longitude_condensed);
end
quiverm(LAND.latitude_condensed(:,:,1),LAND.longitude_condensed(:,:,1),...
    nanmean(LAND.v10_condensed(:,:,n:12:end),3),...
    nanmean(LAND.u10_condensed(:,:,n:12:end),3),'k');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([-2 2]); colormap(cmocean('balance','pivot',0));
c.Label.String = 'FCO_{2} (mol m^{-2} yr^{-1})';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/deltapCO2_RF_Mar.jpg');

% Capture the plot as an image 
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if n == 1 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append'); 
end

end

% Plot month of maximum pCO2 for LAND product
LAND.monthmax = nan(size(LAND.pCO2,1),size(LAND.pCO2,2));
for a = 1:size(LAND.pCO2,1)
    for b = 1:size(LAND.pCO2,2)
        if sum(isnan(squeeze(LAND.pCO2(a,b,:)))) == 12
            LAND.monthmax(a,b) = NaN;
        else
            LAND.monthmax(a,b) = find(squeeze(LAND.pCO2(a,b,:)) == max(squeeze(LAND.pCO2(a,b,:))),1);
        end
    end        
end
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
    LAND.monthmax,...
    0.5:1:12.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([0.5 12.5]); colormap(jet(12));
c.Label.String = 'Month of Maximum pCO_{2}';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/monthmax_LAND.jpg');


%% Plot coastal Laruelle estimates
importLAR

%% Plot random forest Regression estimates

% Plot annual mean RF-predicted pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_RF_clim(:,:,11),3,'omitnan'),...
    300:10:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
%cptcmap('GMT_ocean','flip',true);
colormap(cmocean('haline'));
caxis([300 440]);
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF.jpg');

% Plot RF-predicted pCO2 amplitude
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    log10(SOCATv2020_grid.pco2_RF_amp),...
    1:0.1:2.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
colormap(cmocean('thermal'));
caxis([1 2.5]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_amp.jpg');

% Plot annual mean RF-predicted pCO2 (training data only)
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_RF_validate,3,'omitnan'),...
    280:1:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([280 440]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_validate.jpg');

% Plot annual mean RF-predicted DELTApCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.delpco2_RF_ERA5,3,'omitnan'),...
    -50:1:100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([-50 100]); colormap(cmocean('balance','pivot',0));
c.Label.String = '\DeltapCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/deltapCO2_RF.jpg');

% Plot RF-predicted pCO2 
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '/Users/sharp/Desktop/pCO2_RF.gif';
MonthNames={'January' 'February' 'March' 'April' 'May' 'June' 'July' 'August' 'September' 'October' 'November' 'December'};
for n = 1:12
    worldmap(latlims,lonlims);
    title(MonthNames(n),'fontsize',18);
    set(gcf,'Position',pos)
    setm(gca,'ffacecolor',ocncol);
    setm(gca,'fontsize',fntsz);
    set(gca,'fontsize',fntsz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
        SOCATv2020_grid.pco2_RF_clim(:,:,n),...
        300:5:440,'LineStyle','none');
    geoshow(land, 'FaceColor',lndcol);
    c=colorbar;
    colormap(cmocean('haline'));
    caxis([300 440]);
    c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
    c.Label.String = 'pCO_{2} (\muatm)';
    c.Label.FontSize = fntsz;
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1); 
    end
end

% Plot RF-predicted CO2 flux
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '/Users/sharp/Desktop/CO2_Flux.gif';
MonthNames={'January' 'February' 'March' 'April' 'May' 'June' 'July' 'August' 'September' 'October' 'November' 'December'};
for n = 1:12
    worldmap(latlims,lonlims);
    title(MonthNames(n),'fontsize',18);
    set(gcf,'Position',pos)
    setm(gca,'ffacecolor',ocncol);
    setm(gca,'fontsize',fntsz);
    set(gca,'fontsize',fntsz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
        mean(SOCATv2020_grid.Fco2_RF_ERA5(:,:,n:12:end),3,'omitnan'),...
        -2:0.05:2,'LineStyle','none');
    SOCATv2020_grid.longitude_condensed = repmat([220:1:255]',1,46);
    SOCATv2020_grid.latitude_condensed = repmat(15:1:60,36,1);
    for i=1:264
        SOCATv2020_grid.v10_condensed(:,:,i) = interp2(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
            SOCATv2020_grid.v10(:,:,i),...
            SOCATv2020_grid.latitude_condensed,SOCATv2020_grid.longitude_condensed);
        SOCATv2020_grid.u10_condensed(:,:,i) = interp2(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
            SOCATv2020_grid.u10(:,:,i),...
            SOCATv2020_grid.latitude_condensed,SOCATv2020_grid.longitude_condensed);
    end
    quiverm(SOCATv2020_grid.latitude_condensed(:,:,1),SOCATv2020_grid.longitude_condensed(:,:,1),...
        nanmean(SOCATv2020_grid.v10_condensed(:,:,n:12:end),3),...
        nanmean(SOCATv2020_grid.u10_condensed(:,:,n:12:end),3),'k');
    geoshow(land, 'FaceColor',lndcol);
    c=colorbar; caxis([-2 2]); colormap(cmocean('balance','pivot',0));
    c.Label.String = 'FCO_{2} (mol m^{-2} yr^{-1})';
    c.Label.FontSize = fntsz;

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
end

% Plot RFR values vs SOCAT grid cells (climatological)
figure;
set(gca,'fontsize',fntsz);
Index = ~isnan(SOCATv2020_grid.all.pco2_ave_weighted_clim) & ~isnan(SOCATv2020_grid.pco2_RF_clim);
scatter(SOCATv2020_grid.all.pco2_ave_weighted_clim(Index),...
    SOCATv2020_grid.pco2_RF_clim(Index),...
    10,log10(SOCATv2020_grid.distance_from_shore(Index)),'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCAT gridded monthly mean pCO_{2} (\muatm)','fontsize',16);
xlim([100 1000]);
ylabel('Sharp (RFR estimated) climatological pCO_{2} (\muatm)','fontsize',16);
ylim([100 1000]);
c=colorbar; colormap(cmocean('balance'));
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/RF_vs_SOCAT_individualmonths_scatter.jpg');


% Plot RFR values vs SOCAT grid cells (monthly)
figure;
set(gca,'fontsize',fntsz);
Index = ~isnan(SOCATv2020_grid.all.pco2_ave_weighted) & ~isnan(SOCATv2020_grid.pco2_RF);
scatter(SOCATv2020_grid.all.pco2_ave_weighted(Index),...
    SOCATv2020_grid.pco2_RF(Index),...
    10,log10(SOCATv2020_grid.distance_from_shore(Index)),'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCAT gridded monthly pCO_{2} (\muatm)','fontsize',16);
xlim([100 1000]);
ylabel('Sharp (RFR estimated) monthly pCO_{2} (\muatm)','fontsize',16);
ylim([100 1000]);
c=colorbar; colormap(cmocean('balance'));
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/RF_vs_SOCAT_individualmonths_scatter.jpg');

% Plot wind vs. delta pCO2 RMSE ratio in flux variability
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.RMSE_ratio_RF_ERA5_dtr,...
    0:0.025:2,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([0 2]); colormap(cmocean('balance','pivot',1));
c.Label.String = '|U|^{2} effect / \DeltapCO_{2} effect';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/Fco2_RMSE_ratio_RF.jpg');

% Plot month of maximum pCO2 for RFR product
SOCATv2020_grid.monthmax = nan(size(SOCATv2020_grid.pco2_RF,1),size(SOCATv2020_grid.pco2_RF,2));
for a = 1:size(SOCATv2020_grid.pco2_RF,1)
    for b = 1:size(SOCATv2020_grid.pco2_RF,2)
        if sum(isnan(squeeze(SOCATv2020_grid.pco2_RF_clim(a,b,:)))) == 12
            SOCATv2020_grid.monthmax(a,b) = NaN;
        else
            SOCATv2020_grid.monthmax(a,b) = find(squeeze(SOCATv2020_grid.pco2_RF_clim(a,b,:)) == max(squeeze(SOCATv2020_grid.pco2_RF_clim(a,b,:))));
        end
    end        
end
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.monthmax,...
    0.5:1:12.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([0.5 12.5]); colormap(jet(12));
c.Label.String = 'Month of Maximum pCO_{2}';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/monthmax_RF.jpg');

% Plot wind vs. delta pCO2 correlation
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.U_vs_delpCO2_corr,...
    -1:0.02:1,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(cmocean('balance','pivot',0)); caxis([-1 1]);
c.Label.String = 'Correlation (Wind speed vs. |\DeltapCO_{2}|)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/U_vs_delpCO2_corr_LAND.jpg');

% Plot season of maximum delta pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.seasmax_delpCO2,...
    0.5:4.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(linspecer(4)); caxis([0.5 4.5]);
c.Label.String = 'Season of maximum \DeltapCO_{2}';
c.Label.FontSize = fntsz;
c.Ticks = [1 2 3 4];
c.TickLabels = {'Winter (DJF)' 'Spring (MAM)' 'Summer (JJA)' 'Autumn (SON)'};
exportgraphics(gcf,'/Users/sharp/Desktop/seas_max_delpCO2.jpg');

% Plot season of maximum wind speed
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.seasmax_wind_speed,...
    0.5:4.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(linspecer(4)); caxis([0.5 4.5]);
c.Label.String = 'Season of maximum wind speed';
c.Label.FontSize = fntsz;
c.Ticks = [1 2 3 4];
c.TickLabels = {'Winter (DJF)' 'Spring (MAM)' 'Summer (JJA)' 'Autumn (SON)'};
exportgraphics(gcf,'/Users/sharp/Desktop/seas_max_wind.jpg');

% Plot season of maximum temperature
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.seasmax_SST,...
    0.5:4.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(linspecer(4)); caxis([0.5 4.5]);
c.Label.String = 'Season of maximum SST';
c.Label.FontSize = fntsz;
c.Ticks = [1 2 3 4];
c.TickLabels = {'Winter (DJF)' 'Spring (MAM)' 'Summer (JJA)' 'Autumn (SON)'};
exportgraphics(gcf,'/Users/sharp/Desktop/seas_max_SST.jpg');

% Plot season of maximum chlorophyll
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.seasmax_CHL,...
    0.5:4.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(linspecer(4)); caxis([0.5 4.5]);
c.Label.String = 'Season of maximum chlorophyll';
c.Label.FontSize = fntsz;
c.Ticks = [1 2 3 4];
c.TickLabels = {'Winter (DJF)' 'Spring (MAM)' 'Summer (JJA)' 'Autumn (SON)'};
exportgraphics(gcf,'/Users/sharp/Desktop/seas_max_CHL.jpg');

%% Plot difference maps

% Plot Landschutzer climatological pCO2 - RF-predicted climatological pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_RF_clim_2015,3,'omitnan')-mean(LAND.pCO2,3,'omitnan'),...
    -50:1:150,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-50 150]);
colormap(cmocean('balance','pivot',0));
c.TickLabels = {'-50' '0' '50' '100' '150+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_minus_LAND.jpg');

% Plot Landschutzer amplitude - RF-predicted amplitude
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.pco2_RF_amp_2015(:,:,1)-LAND.pCO2_amp(:,:,1),...
    -50:1:150,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-50 150]);
colormap(cmocean('balance','pivot',0));
c.TickLabels = {'-50' '0' '50' '100' '150+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_minus_LAND_Amp.jpg');

% Plot Laruelle climatological pCO2 - RF-predicted climatological pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_RF_clim_2015,3,'omitnan')-mean(LAR.pCO2,3,'omitnan'),...
    -50:1:150,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-50 150]);
colormap(cmocean('balance','pivot',0));
c.TickLabels = {'-50' '0' '50' '100' '150+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_minus_LAR.jpg');

% Plot Laruelle amplitude - RF-predicted amplitude
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.pco2_RF_amp_2015(:,:,1)-LAR.pCO2_amp(:,:,1),...
    -50:1:150,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-50 150]);
colormap(cmocean('balance','pivot',0));
c.TickLabels = {'-50' '0' '50' '100' '150+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF_minus_LAR_Amp.jpg');

%% Plot SOM-FFN estimates

% Plot annual mean FFN-predicted pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_FFN,3,'omitnan'),...
    280:1:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([280 440]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_FFN.jpg');

% Plot annual mean FFN-predicted pCO2 (training data only)
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_FFN_validate,3,'omitnan'),...
    280:1:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([280 440]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_FFN_validate.jpg');

%% Plots of predictor variables

% Plot SST
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    nanmean(SOCATv2020_grid.SST,3),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet);
c.Label.String = 'Temperature (degC)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/SST.jpg');

% Plot SSS
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    nanmean(SOCATv2020_grid.SSS,3),...
    28:0.125:36,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([28 36]);
c.Label.String = 'Salinity';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/SSS.jpg');

% Plot SSH
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    nanmean(SOCATv2020_grid.SSH,3),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); %caxis([28 36]);
c.Label.String = 'Sea Surface Height';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/SSH.jpg');

% Plot CHL
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    log10(nanmean(SOCATv2020_grid.CHL(:,:,12*20),3)),...
    -2:0.05:2,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([-2 2]);
c.Label.String = 'Chlorophyll-a (mg m^{-3})';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/CHL.jpg');

% Plot Wind Speed
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    nanmean(SOCATv2020_grid.wind_speed,3),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet);
c.Label.String = 'Wind Speed (m/s)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/WindSpeed.jpg');

% Plot Mixed Layer Depth
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    nanmean(SOCATv2020_grid.MLD,3),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet);
c.Label.String = 'Mixed Layer Depth (m)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/MLD.jpg');

% Plot Bottom Depth
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    SOCATv2020_grid.bottomdepth(:,:,1),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(jet); caxis([0 6000]);
c.Label.String = 'Bottom Depth (m)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/Bathymetry.jpg');

%% Plot mooring-related figures

% Landshutzer/mooring amplitude map
figure
worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),LAND.pCO2_amp,0:1:100,'LineStyle','none');
geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
for n=1:numel(moornames)
    % Find closest point in Landschutzer climatology to lat-lon
    latidx = ...
        find(abs(LAND.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(LAND.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(LAND.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(LAND.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    %if strcmp('Southeast',moornames{n}); latidx = latidx - 1; end
    scatterm(LAND.latitude(lonidx,latidx,1),LAND.longitude(lonidx,latidx,1),...
        200,MOORING.(moornames{n}).pCO2SW_amplitude,'filled',...
        'MarkerEdgeColor','k','LineWidth',2);
    if n==7
    textm(LAND.latitude(lonidx,latidx,1)+1,LAND.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','white','FontWeight','bold','fontsize',16);
    else
    textm(LAND.latitude(lonidx,latidx,1),LAND.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','white','FontWeight','bold','fontsize',16);
    end
end
c=colorbar; colormap(jet); caxis([0 100]);
c.Label.String = 'pCO_{2} (\muatm)';
c.FontSize = 18;
exportgraphics(gcf,'/Users/sharp/Desktop/Land_mooring_amp_map.jpg');

% Landshutzer/mooring amplitude map
figure
worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAND.latitude(:,:,1),LAND.longitude(:,:,1),...
    log10(LAND.pCO2_amp),1:0.1:2.5,'linestyle','none');
geoshow(land, 'FaceColor',lndcol);
for n=1:numel(moornames)
    % Find closest point in RF climatology to lat-lon
    latidx = ...
        find(abs(LAND.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(LAND.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(LAND.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(LAND.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    scatterm(LAND.latitude(lonidx,latidx,1),LAND.longitude(lonidx,latidx,1),...
        400,log10(MOORING.(moornames{n}).pCO2SW_amplitude),'filled',...
        'MarkerEdgeColor','k','LineWidth',4);
    if n==6
    textm(LAND.latitude(lonidx,latidx,1)+1,LAND.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    elseif n==1
    textm(LAND.latitude(lonidx,latidx,1),LAND.longitude(lonidx,latidx,1)+6,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    elseif n==2
    textm(LAND.latitude(lonidx,latidx,1)+1,LAND.longitude(lonidx,latidx,1)+2,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    else
    textm(LAND.latitude(lonidx,latidx,1),LAND.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    end
end
c=colorbar;
colormap(cmocean('thermal'));
caxis([1 2.5]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/Land_mooring_amp_map.jpg');

% Laruelle/mooring amplitude map
figure
worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(LAR.latitude(:,:,1),LAR.longitude(:,:,1),LAR.pCO2_amp,...
    0:1:100,'LineStyle','none');
geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
for n=1:numel(moornames)
    % Find closest point in Laruelle climatology to lat-lon
    latidx = ...
        find(abs(LAR.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(LAR.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(LAR.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(LAR.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    %if strcmp('Southeast',moornames{n}); latidx = latidx - 1; end
    scatterm(LAR.latitude(lonidx,latidx,1),LAR.longitude(lonidx,latidx,1),...
        200,MOORING.(moornames{n}).pCO2SW_amplitude,'filled',...
        'MarkerEdgeColor','k','LineWidth',2);
    if n==7
    textm(LAR.latitude(lonidx,latidx,1)+1,LAR.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','white','FontWeight','bold','fontsize',16);
    else
    textm(LAR.latitude(lonidx,latidx,1),LAR.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','white','FontWeight','bold','fontsize',16);
    end
end
c=colorbar; colormap(jet); caxis([0 100]);
c.Label.String = 'pCO_{2} (\muatm)';
c.FontSize = 18;
exportgraphics(gcf,'/Users/sharp/Desktop/Lar_mooring_amp_map.jpg');


% RFR/mooring amplitude map
figure
worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    log10(SOCATv2020_grid.pco2_RF_amp),1:0.1:2.5,'linestyle','none');
geoshow(land, 'FaceColor',lndcol);
for n=1:numel(moornames)
    % Find closest point in RF climatology to lat-lon
    latidx = ...
        find(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    scatterm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1),...
        400,log10(MOORING.(moornames{n}).pCO2SW_amplitude),'filled',...
        'MarkerEdgeColor','k','LineWidth',4);
    if n==6
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1)+1,SOCATv2020_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    elseif n==1
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1)+6,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    elseif n==2
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1)+1,SOCATv2020_grid.longitude(lonidx,latidx,1)+2,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    else
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    end
end
c=colorbar;
colormap(cmocean('thermal'));
caxis([1 2.5]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/RF_mooring_amp_map.jpg');

% RFR/mooring mean map
figure
worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pco2_RF,3,'omitnan'),300:10:440,'linestyle','none');
geoshow(land, 'FaceColor',lndcol);
for n=1:numel(moornames)
    % Find closest point in RF climatology to lat-lon
    latidx = ...
        find(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    scatterm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1),...
        400,mean(MOORING.(moornames{n}).pCO2SW_monthly(:,1),1),'filled',...
        'MarkerEdgeColor','k','LineWidth',4);
    if n==7
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1)+1,SOCATv2020_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    elseif n==1
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1)+6,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    elseif n==2
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1)+1,SOCATv2020_grid.longitude(lonidx,latidx,1)+2,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
    else
    textm(SOCATv2020_grid.latitude(lonidx,latidx,1),SOCATv2020_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
    end
end
c=colorbar;
colormap(cmocean('haline'));
caxis([300 440]);
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/RF_mooring_mean_validate_map.jpg');



% Plot climatology pCO2
figure
set(gcf,'Position',[100, 400, 1680, 420])
set(gca,'fontsize',16)
hold on
date = datenum([repmat(1998,264,1) SOCATv2020_grid.month_since_1998 ones(264,1)]);
p1=plot(date,squeeze(mean(mean(SOCATv2020_grid.pco2_RF_validate,1,'omitnan'),2,'omitnan')),'linewidth',2);
p2=plot(date,squeeze(mean(mean(SOCATv2020_grid.pco2_RF,1,'omitnan'),2,'omitnan')),'linewidth',2);
%p3=plot(date(1:216),squeeze(mean(mean(LAR.pCO2,1,'omitnan'),2,'omitnan')),'linewidth',2);
legend([p1 p2],{'Sharp, RFR (training data)'...
                   'Sharp, RFR (all data)'},...
                   'location','southeast');
ylabel('pCO_{2} (\muatm)','fontsize',20);
xlabel('Year','fontsize',20);
datetick('x','yyyy');


%% Plot TA
% Plot annual mean RF-predicted pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.TA_ESPER_LIR(:,:,:),3,'omitnan'),...
    100,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
%cptcmap('GMT_ocean','flip',true);
colormap(cmocean('haline'));
%caxis([300 440]);
%c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'TA (\mumol kg^{-1})';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/TA.jpg');

%% Plot pH
% Plot annual mean pH from RF-predicted pCO2 and ESPER-predicted TA
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.pH_clim,3,'omitnan'),...
    8:0.001:8.1,'LineStyle','none');
% [C,h]=contourm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
%     mean(SOCATv2020_grid.pH_clim,3,'omitnan'),...
%     8:.01:8.1,'LineStyle','-','linecolor','k','linewidth',2);
% clabelm(C,h,[8:0.05:8.1]);
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
colormap(cmocean('haline'));
caxis([8 8.1]);
%c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'pH';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/pH.jpg');

%% Plot Omega
% Plot annual mean Omega Aragonite from RF-predicted pCO2 and ESPER-predicted TA
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.OmA_clim,3,'omitnan'),...
    1:.025:4,'LineStyle','none');
[C,h]=contourm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.OmA_clim,3,'omitnan'),...
    1:.5:4,'LineStyle','-','linecolor','k','linewidth',2);
clabelm(C,h,[1 1.5 2 2.5 3 3.5 4]);
% scatterm(SOCATv2020_grid.latitude(lonidx,latidx,1),...
%     SOCATv2020_grid.longitude(lonidx,latidx,1),200,'^k','filled');
% scatterm(SOCATv2020_grid.latitude(lonidx2,latidx2,1),...
%     SOCATv2020_grid.longitude(lonidx2,latidx2,1),200,'^k','filled');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
%cptcmap('GMT_ocean','flip',true);
colormap(cmocean('deep'));
caxis([1 4]);
%c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = '\Omega_{A}';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/OmA.jpg');
