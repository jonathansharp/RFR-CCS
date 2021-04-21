% Figure 2
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol  = [0.94 0.97 1.0];
fntsz   = 12;
lndcol  = [0.2 0.2 0.2];

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
        mean(SOCATv2020_grid.pco2_RF(:,:,n:12:end),3,'omitnan'),...
        300:10:440,'LineStyle','none');
    geoshow(land, 'FaceColor',lndcol);
    colormap(cmocean('haline'));
    caxis([300 440]);
    
end

% Define axis
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',20);
c=colorbar;
caxis([300 440]);
colormap(cmocean('haline',14));
c.Ticks = [300 320 340 360 380 400 420 440];
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = 20;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/Figure2.jpg');
