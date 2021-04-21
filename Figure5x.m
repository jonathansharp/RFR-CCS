% Figure 5
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol = [0.94 0.97 1.0];
fntsz = 22;
lndcol = [0.2 0.2 0.2];

% Import SOCATv4 gridded files
% importSOCATv4gridded
% importSOCATv4gridded_coastal
SOCATv4_gridded.distance_from_shore = ...
    dist2coast(SOCATv4_gridded.lat,SOCATv4_gridded.lon);
SOCATv4_gridded.distance_from_shore = ...
    repmat(SOCATv4_gridded.distance_from_shore,1,1,...
    size(SOCATv4_gridded.month_since_1998,1));

% Initialize figure
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
tl = axes('Position',[0.05 0.56 0.4 0.425],'Box','on');
tr = axes('Position',[0.55 0.56 0.4 0.425],'Box','on');
bl = axes('Position',[0.05 0.06 0.4 0.425],'Box','on');
br = axes('Position',[0.55 0.06 0.4 0.425],'Box','on');

% Interpolate values to open-ocean grid cells
for m = 1:12
    % RFR climatology
    interp = griddedInterpolant(SOCATv2020_grid.lon,SOCATv2020_grid.lat,SOCATv2020_grid.pco2_RF_clim_2015(:,:,m));
    SOCATv2020_grid.pco2_RF_clim_open_2015(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    % Landschutzer product
    interp = griddedInterpolant(SOCATv2020_grid.lon,SOCATv2020_grid.lat,LAND.pCO2(:,:,m));
    LAND.pCO2_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    
end
for m = 1:264
    % RFR monthly product
    interp = griddedInterpolant(SOCATv2020_grid.lon,SOCATv2020_grid.lat,SOCATv2020_grid.pco2_RF(:,:,m));
    SOCATv2020_grid.pco2_RF_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
end

% Plot RFR values vs SOCAT grid cells (climatological)
axes(tl);
set(gca,'fontsize',fntsz);
IndexOpen = ~isnan(SOCATv4_gridded.pCO2_mon_mean_clim) & ~isnan(SOCATv2020_grid.pco2_RF_clim_open_2015);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean_clim) & ~isnan(SOCATv2020_grid.pco2_RF_clim_2015);
scatter([SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen);SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast)],...
    [SOCATv2020_grid.pco2_RF_clim_open_2015(IndexOpen);SOCATv2020_grid.pco2_RF_clim_2015(IndexCoast)],5,...
    [log10(SOCATv4_gridded.distance_from_shore(IndexOpen));log10(SOCATv2020_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean pCO_{2} (\muatm)','fontsize',fntsz);
xlim([100 900]);
ylabel('RFR-NEP-clim pCO_{2} (\muatm)','fontsize',fntsz);
ylim([100 800]);
c=colorbar;
colormap(cmocean('balance'));
caxis([1 3]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast);SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen)],...
    [SOCATv2020_grid.pco2_RF_clim_2015(IndexCoast);SOCATv2020_grid.pco2_RF_clim_open_2015(IndexOpen)]),2))),'fontsize',fntsz);
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;

% Plot Landschutzer values vs SOCAT grid cells
axes(tr);
set(gca,'fontsize',fntsz);
IndexOpen = ~isnan(SOCATv4_gridded.pCO2_mon_mean_clim) & ~isnan(LAND.pCO2_open);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean_clim) & ~isnan(LAND.pCO2);
scatter([SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen);SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast)],...
    [LAND.pCO2_open(IndexOpen);LAND.pCO2(IndexCoast)],5,...
    [log10(SOCATv4_gridded.distance_from_shore(IndexOpen));log10(SOCATv2020_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean pCO_{2} (\muatm)','fontsize',fntsz);
xlim([100 900]);
ylabel('Landschützer et al. (2020) pCO_{2} (\muatm)','fontsize',fntsz);
ylim([100 800]);
c=colorbar;
colormap(cmocean('balance'));
caxis([1 3]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast);SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen)],...
    [LAND.pCO2(IndexCoast);LAND.pCO2_open(IndexOpen)]),2))),'fontsize',fntsz);
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;

% Plot RFR values vs SOCAT grid cells (monthly)
axes(bl);
set(gca,'fontsize',fntsz);
IndexOpen = ~isnan(SOCATv4_gridded.pCO2_mon_mean) & ~isnan(SOCATv2020_grid.pco2_RF_open(:,:,1:216));
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean) & ~isnan(SOCATv2020_grid.pco2_RF(:,:,1:216));
scatter([SOCATv4_gridded.pCO2_mon_mean(IndexOpen);SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast)],...
    [SOCATv2020_grid.pco2_RF_open(IndexOpen);SOCATv2020_grid.pco2_RF(IndexCoast)],5,...
    [log10(SOCATv4_gridded.distance_from_shore(IndexOpen));log10(SOCATv2020_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly pCO_{2} (\muatm)','fontsize',fntsz);
xlim([100 900]);
ylabel('RFR-NEP pCO_{2} (\muatm)','fontsize',fntsz);
ylim([100 800]);
c=colorbar;
colormap(cmocean('balance'));
caxis([1 3]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast);SOCATv4_gridded.pCO2_mon_mean(IndexOpen)],...
    [SOCATv2020_grid.pco2_RF(IndexCoast);SOCATv2020_grid.pco2_RF_open(IndexOpen)]),2))),'fontsize',fntsz);
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;

% Plot Laruelle values vs SOCAT grid cells
axes(br);
set(gca,'fontsize',fntsz);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean) & ~isnan(LAR.pCO2);
scatter(SOCATv4_gridded_coastal.pCO2_mon_mean(Index),...
    LAR.pCO2(Index),5,log10(SOCATv2020_grid.distance_from_shore(Index)),'filled','MarkerEdgeColor','none'); hold on;
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly pCO_{2} (\muatm)','fontsize',fntsz);
xlim([100 900]);
ylabel('Laruelle et al. (2017) pCO_{2} (\muatm)','fontsize',fntsz);
ylim([100 800]);
c=colorbar;
colormap(cmocean('balance'));
caxis([1 3]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr(SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast),...
    LAR.pCO2(IndexCoast)),2))),'fontsize',fntsz);
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = fntsz;

exportgraphics(gcf,'/Users/sharp/Desktop/Figure5.jpg');