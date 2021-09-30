%% This script produces Figure A2 from Sharp et al. (in prep)

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
    mean(SOCATv2021_grid.CHL,3,'omitnan'));
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([0 2]);
colormap(parula(8));
c.Ticks = [0 0.5 1 1.5 2];
c.TickLabels = {'0.0' '0.5' '1.0' '1.5' '2.0+'};
c.Label.String = ['Chlorophyll-a (mg m^{-2})'];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'Figures/FigureA2.png');