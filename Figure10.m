%% This script produces Figure 10 from Sharp et al. (in prep)
% Plots relative effests of delta pCO2 and wind speed on CO2 flux
% variability.

%% Initialize figure
figure;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
l = axes('Position',[0.05 0.05 0.4 0.9],'Box','on');
r = axes('Position',[0.55 0.05 0.4 0.9],'Box','on');

ocncol  = [1 1 1];
lndcol  = [1 1 1];
fntsz   = 18;

%% Plot delta pCO2 effect
axes(l); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    pco2_effect,-2:0.2:2,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-2 2]);
colormap(cmocean('balance',20,'pivot',0));
c.Label.String = '{\it\beta}_{pCO2}';
c.Label.FontSize = 28;

%% Plot wind effect
axes(r); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    wind_effect,-2:0.2:2,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-2 2]);
colormap(cmocean('balance',20,'pivot',0));
c.Label.String = '{\it\beta}_{wind}';
c.Label.FontSize = 28;

%% Add figure labels
annotation('textbox',[.03 .85 .1 .1],'String','a','EdgeColor','none','fontsize',32)
annotation('textbox',[.53 .85 .1 .1],'String','b','EdgeColor','none','fontsize',32)

%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure10.jpg');
