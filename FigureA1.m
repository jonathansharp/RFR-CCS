%% This script produces Figure A1 from Sharp et al. (in prep)

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol = [1 1 1];
lndcol = [1 1 1];
fntsz = 18;

% Plot annual mean gridded temp
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.SST,3,'omitnan'));
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([5 30]);
colormap(parula(10));
%c.TickLabels = {'5'};
c.Label.String = ['SST (' char(176) 'C)'];
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/FigureA1.jpg');