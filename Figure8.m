%% This script produces Figure 8 from Sharp et al. (in prep)
% It gives monthly mean fields of CO2 flux determined by from the RFR-NEP
% surface pCO2 product

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol  = [1 1 1];
fntsz   = 12;
lndcol  = [1 1 1];

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 0.7]);

mnth = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
x_ax = [0.02 0.17 0.32 0.47 0.62 0.77 0.02 0.17 0.32 0.47 0.62 0.77];
y_ax = [0.55 0.55 0.55 0.55 0.55 0.55 0.05 0.05 0.05 0.05 0.05 0.05];

for n = 1:12

    % Define axis
    axes('Position',[x_ax(n) y_ax(n) 0.15 0.4],'Box','off');
    
    % Plot monthly RF-predicted pCO2
    worldmap(latlims,lonlims);
    title(mnth{n});
    setm(gca,'ffacecolor',ocncol);
    setm(gca,'fontsize',fntsz);
    set(gca,'fontsize',fntsz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
        mean(SOCATv2020_grid.Fco2_RF_ERA5(:,:,n:12:end),3,'omitnan').*12,...
        -5:0.5:5,'LineStyle','none');
    % Get condensed lat-lon matrices
    SOCATv2020_grid.longitude_condensed = repmat([220:2.5:255]',1,19);
    SOCATv2020_grid.latitude_condensed = repmat(15:2.5:60,15,1);
    % interpolate wind speeds onto condensed grid
    SOCATv2020_grid.v10_condensed = nan(15,19,264);
    SOCATv2020_grid.u10_condensed = nan(15,19,264);
    for i=1:264
        SOCATv2020_grid.v10_condensed(:,:,i) = interp2(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
            SOCATv2020_grid.v10(:,:,i),...
            SOCATv2020_grid.latitude_condensed,SOCATv2020_grid.longitude_condensed);
        SOCATv2020_grid.u10_condensed(:,:,i) = interp2(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
            SOCATv2020_grid.u10(:,:,i),...
            SOCATv2020_grid.latitude_condensed,SOCATv2020_grid.longitude_condensed);
    end
    % Plot wind speed vectors
    quiverm(SOCATv2020_grid.latitude_condensed(:,:,1),SOCATv2020_grid.longitude_condensed(:,:,1),...
        nanmean(SOCATv2020_grid.v10_condensed(:,:,n:12:end),3),...
        nanmean(SOCATv2020_grid.u10_condensed(:,:,n:12:end),3),'k');
    
    geoshow(land, 'FaceColor',lndcol);
    caxis([-5 5]);
    colormap(cmocean('balance',20,'pivot',0));
    
end

% Define axis
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',20);
c=colorbar;
caxis([-5 5]);
colormap(cmocean('balance',20,'pivot',0));
c.Ticks = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
c.TickLabels = {'-5' '-4' '-3' '-2' '-1' '0' '1' '2' '3' '4' '5'};
c.Label.String = '{\itF}_{CO2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 20;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/Figure8.jpg');
