%% This script produces Figure A4 from Sharp et al. (in prep)

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol = [1 1 1];
lndcol = [0.5 0.5 0.5];
fntsz = 18;

% Plot annual mean gridded temp
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.wind_speed,3,'omitnan'));
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([0 10]);
colormap(parula(10));
c.Label.String = ['Wind Speed at 10m (m s^{-1})'];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'Figures/FigureA4.png');