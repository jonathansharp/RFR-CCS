%% This script produces Figure 4 from Sharp et al. (in prep)

% Import SOCATv4 gridded files
importSOCATv4gridded
importSOCATv4gridded_coastal

SOCATv4_gridded.distance_from_shore = ...
    dist2coast(SOCATv4_gridded.lat,SOCATv4_gridded.lon);
SOCATv4_gridded.distance_from_shore = ...
    repmat(SOCATv4_gridded.distance_from_shore,1,1,...
    size(SOCATv4_gridded.month_since_1998,1));

%% Initialize figure

labelsz = 24;
fontsz = 20;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.75 0.75])
tl = axes('Position',[0.06 0.59 0.37 0.38],'Box','on');
tr = axes('Position',[0.53 0.59 0.37 0.38],'Box','on');
bl = axes('Position',[0.06 0.09 0.37 0.38],'Box','on');
br = axes('Position',[0.53 0.09 0.37 0.38],'Box','on');

%% Interpolate values to open-ocean grid cells
for m = 1:12
    % RFR climatology
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,SOCATv2021_grid.pco2_RF_clim_2015(:,:,m));
    SOCATv2021_grid.pco2_RF_clim_open_2015(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    % Landschutzer product
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,LAND.pCO2(:,:,m));
    LAND.pCO2_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    
end
for m = 1:264
    % RFR monthly product
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,SOCATv2021_grid.pco2_RF(:,:,m));
    SOCATv2021_grid.pco2_RF_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
end

%% Plot RFR values vs SOCAT grid cells (climatological)
axes(tl);
IndexOpen = ~isnan(SOCATv4_gridded.pCO2_mon_mean_clim) & ~isnan(SOCATv2021_grid.pco2_RF_clim_open_2015);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean_clim) & ~isnan(SOCATv2021_grid.pco2_RF_clim_2015);
scatter([SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen);SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast)],...
    [SOCATv2021_grid.pco2_RF_clim_open_2015(IndexOpen);SOCATv2021_grid.pco2_RF_clim_2015(IndexCoast)],20,...
    [log10(SOCATv4_gridded.distance_from_shore(IndexOpen));log10(SOCATv2021_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([100 900]);
ylabel('RFR-CCS-clim {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([100 800]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
xticks([200 400 600 800]);
yticks([200 400 600 800]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast);SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen)],...
    [SOCATv2021_grid.pco2_RF_clim_2015(IndexCoast);SOCATv2021_grid.pco2_RF_clim_open_2015(IndexOpen)]),2))),'fontsize',labelsz);

sum(~isnan([SOCATv2021_grid.pco2_RF_clim_open_2015(IndexOpen);SOCATv2021_grid.pco2_RF_clim_2015(IndexCoast)]))

%% Plot Landschutzer values vs SOCAT grid cells
axes(tr);
IndexOpen = ~isnan(SOCATv4_gridded.pCO2_mon_mean_clim) & ~isnan(LAND.pCO2_open);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean_clim) & ~isnan(LAND.pCO2);
scatter([SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen);SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast)],...
    [LAND.pCO2_open(IndexOpen);LAND.pCO2(IndexCoast)],20,...
    [log10(SOCATv4_gridded.distance_from_shore(IndexOpen));log10(SOCATv2021_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([100 900]);
ylabel('L20 {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([100 800]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
xticks([200 400 600 800]);
yticks([200 400 600 800]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean_clim(IndexCoast);SOCATv4_gridded.pCO2_mon_mean_clim(IndexOpen)],...
    [LAND.pCO2(IndexCoast);LAND.pCO2_open(IndexOpen)]),2))),'fontsize',labelsz);

sum(~isnan([LAND.pCO2_open(IndexOpen);LAND.pCO2(IndexCoast)]))

%% Plot RFR values vs SOCAT grid cells (monthly coastal)
axes(bl);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean) & ~isnan(SOCATv2021_grid.pco2_RF(:,:,1:216)) & ~isnan(LAR.pCO2);
scatter([SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast)],...
    [SOCATv2021_grid.pco2_RF(IndexCoast)],20,...
    [log10(SOCATv2021_grid.distance_from_shore(IndexCoast))],...
    'filled','MarkerEdgeColor','none'); hold on;
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([100 900]);
ylabel('RFR-CCS {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([100 800]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
xticks([200 400 600 800]);
yticks([200 400 600 800]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr([SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast)],...
    [SOCATv2021_grid.pco2_RF(IndexCoast)]),2))),'fontsize',labelsz);

sum(~isnan([SOCATv2021_grid.pco2_RF(IndexCoast)]))

%% Plot Laruelle values vs SOCAT grid cells
axes(br);
IndexCoast = ~isnan(SOCATv4_gridded_coastal.pCO2_mon_mean) & ~isnan(LAR.pCO2);
scatter(SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast),...
    LAR.pCO2(IndexCoast),20,log10(SOCATv2021_grid.distance_from_shore(IndexCoast)),'filled','MarkerEdgeColor','none'); hold on;
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([100 900]);
ylabel('L17 {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([100 800]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
xticks([200 400 600 800]);
yticks([200 400 600 800]);
text(200,700,strcat('R^{2} =',{' '},...
    num2str(round(corr(SOCATv4_gridded_coastal.pCO2_mon_mean(IndexCoast),...
    LAR.pCO2(IndexCoast)),2))),'fontsize',labelsz);

sum(~isnan(LAR.pCO2(IndexCoast)))

%% Colorbar
ax=axes('Position',[0.06 0.05 0.9 0.9],'Box','off');
set(gca,'fontsize',fontsz);
c=colorbar;
caxis([1 3]);
colormap(cmocean('haline'));
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = labelsz;
ax.Visible = 'off';

annotation('textbox',[0.4, 0.95, 0, 0],'string','a','fontsize',24)
annotation('textbox',[0.86, 0.95, 0, 0],'string','b','fontsize',24)
annotation('textbox',[0.4, 0.45, 0, 0],'string','c','fontsize',24)
annotation('textbox',[0.86, 0.45, 0, 0],'string','d','fontsize',24)

exportgraphics(gcf,'Figures/Figure4.png');