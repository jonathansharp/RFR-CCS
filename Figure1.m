%% This script produces Figure 1 from Sharp et al. (in prep)

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol = [1 1 1];
lndcol = [1 1 1];
fntsz = 18;

% Plot annual mean gridded pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.all.pco2_ave_weighted_clim,3,'omitnan'));
geoshow(land, 'FaceColor',lndcol,'linewidth',1);
c=colorbar;
caxis([300 440]);
colormap(cmocean('haline',14));
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = '{\itp}CO_{2(sw)} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/Figure1.jpg');