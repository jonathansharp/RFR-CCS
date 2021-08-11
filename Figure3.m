%% This script produces Figure 3 from Sharp et al. (in prep)

%% Initialize figure
figure;
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
tl = axes('Position',[0.02 0.5 0.425 0.4],'Box','on');
tr = axes('Position',[0.42 0.5 0.425 0.4],'Box','on');
bl = axes('Position',[0.02 0.05 0.425 0.4],'Box','on');
br = axes('Position',[0.42 0.05 0.425 0.4],'Box','on');

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol = [1 1 1];
lndcol = [0.5 0.5 0.5];
fntsz = 18;

%% Plot Landschutzer climatological pCO2 - RF-predicted climatological pCO2 (Xa)
axes(tl); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    mean(SOCATv2021_grid.pco2_RF_clim_2015,3,'omitnan')-mean(LAND.pCO2,3,'omitnan'));
geoshow(land,'FaceColor',lndcol,'linestyle','none');
caxis([-50 50]);
textm(50,245,'L20','fontsize',24,'color','w','fontweight','bold');
colormap(cmocean('balance','pivot',0));

%% Plot Landschutzer amplitude - RF-predicted amplitude (Xb)
axes(tr); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    SOCATv2021_grid.pco2_RF_amp_2015-LAND.pCO2_amp);
geoshow(land,'FaceColor',lndcol,'linestyle','none');
caxis([-50 50]);
textm(50,245,'L20','fontsize',24,'color','w','fontweight','bold');
colormap(cmocean('balance','pivot',0));

%% Plot Laruelle climatological pCO2 - RF-predicted climatological pCO2
axes(bl); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    mean(SOCATv2021_grid.pco2_RF_clim_2015,3,'omitnan')-mean(LAR.pCO2,3,'omitnan'));
geoshow(land,'FaceColor',lndcol,'linestyle','none');
caxis([-50 50]);
textm(50,245,'L17','fontsize',24,'color','w','fontweight','bold');
colormap(cmocean('balance','pivot',0));

%% Plot Laruelle amplitude - RF-predicted amplitude
axes(br); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    SOCATv2021_grid.pco2_RF_amp_2015-LAR.pCO2_amp);
geoshow(land,'FaceColor',lndcol,'linestyle','none');
caxis([-50 50]);
textm(50,245,'L17','fontsize',24,'color','w','fontweight','bold');
colormap(cmocean('balance','pivot',0));

%% Colorbar
ax=axes('Position',[0 0.05 0.9 0.9],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
caxis([-50 50]);
colormap(cmocean('balance','pivot',0));
c.Ticks = [-50 -25 0 25 50];
c.TickLabels = {'-50' '-25' '0' '25' '50'};
c.Label.String = '\Delta{\itp}CO_{2(sw)} (\muatm)';
c.Label.FontSize = 24;
ax.Visible = 'off';

annotation('textbox',[0.02, 0.92, 0, 0],'string','a','fontsize',24)
annotation('textbox',[0.42, 0.92, 0, 0],'string','b','fontsize',24)
annotation('textbox',[0.02, 0.445, 0, 0],'string','c','fontsize',24)
annotation('textbox',[0.42, 0.445, 0, 0],'string','d','fontsize',24)

annotation('textbox',[0.07, 0.98, 0.4, 0],'string','Annual mean differences',...
    'fontsize',24,'EdgeColor','none');
annotation('textbox',[0.44, 0.98, 0.4, 0],'string','Seasonal amplitude differences',...
    'fontsize',24,'EdgeColor','none')

%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure3.jpg');

%% Area-weighted differences between RFR-NEP vs. L17 and L20
% area weights
area_weights_open = SOCATv2021_grid.area_km2;
area_weights_open(isnan(LAND.pCO2(:,:,1))) = NaN;
area_weights_coast = SOCATv2021_grid.area_km2;
area_weights_coast(isnan(LAR.pCO2(:,:,1))) = NaN;

% Mean: RFR-NEP - L20
RFR_NEP_Clim_L20 = ...
    sum(sum((mean(SOCATv2021_grid.pco2_RF_clim_2015,3,'omitnan')-...
            mean(LAND.pCO2,3,'omitnan')).*...
            area_weights_open,'omitnan'),'omitnan')./...
            sum(sum(area_weights_open,'omitnan'),'omitnan')

% Amplitude: RFR-NEP - L20
RFR_NEP_Clim_L20_amp = ...
    sum(sum((SOCATv2021_grid.pco2_RF_amp_2015-LAND.pCO2_amp).*...
            area_weights_open,'omitnan'),'omitnan')./...
            sum(sum(area_weights_open,'omitnan'),'omitnan')

% Mean: RFR-NEP - L17
RFR_NEP_Clim_L17 = ...
    sum(sum((mean(SOCATv2021_grid.pco2_RF_clim_2015,3,'omitnan')-...
            mean(LAR.pCO2,3,'omitnan')).*...
            area_weights_coast,'omitnan'),'omitnan')./...
            sum(sum(area_weights_coast,'omitnan'),'omitnan')

% Amplitude: RFR-NEP - L17
RFR_NEP_Clim_L17_amp = ...
    sum(sum((SOCATv2021_grid.pco2_RF_amp_2015-LAR.pCO2_amp).*...
            area_weights_coast,'omitnan'),'omitnan')./...
            sum(sum(area_weights_coast,'omitnan'),'omitnan')

