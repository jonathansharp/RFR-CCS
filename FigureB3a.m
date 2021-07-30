%% This script produces Figure SX from Sharp et al. (in prep)
% 

ocncol  = [1 1 1];
lndcol  = [1 1 1];
fntsz   = 12;

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
    pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
        SOCATv2021_grid.pco2_RF_clim_2015(:,:,n)-LAND.pCO2(:,:,n),...
        'LineStyle','none');
    geoshow(land, 'FaceColor',lndcol);
    caxis([-50 50]);
    colormap(cmocean('balance','pivot',0));
    
end

%% Colorbar
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',20);
c=colorbar;
caxis([-50 50]);
colormap(cmocean('balance','pivot',0));
c.Ticks = [-50 -25 0 25 50];
c.TickLabels = {'-50' '-25' '0' '25' '50+'};
c.Label.String = '\Delta{\itp}CO_{2(sw)} (\muatm)';
c.Label.FontSize = 20;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/FigureB3a.jpg');
