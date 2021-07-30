%% Figure 13

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol  = [0.94 0.97 1.0];
fntsz   = 12;
lndcol  = [0.2 0.2 0.2];
cmap    = jet;

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 0.7]);

mnth = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
x_ax = [0.02 0.17 0.32 0.47 0.62 0.77 0.02 0.17 0.32 0.47 0.62 0.77];
y_ax = [0.55 0.55 0.55 0.55 0.55 0.55 0.05 0.05 0.05 0.05 0.05 0.05];

for n = 1:12

    % Define axis
    axes('Position',[x_ax(n) y_ax(n) 0.15 0.4],'Box','off');
    
    % Plot interannual variability in FCO2
    worldmap(latlims,lonlims);
    title(mnth{n});
    setm(gca,'ffacecolor',ocncol);
    setm(gca,'fontsize',fntsz);
    set(gca,'fontsize',fntsz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
        SOCATv2020_grid.Fco2_monthly_std(:,:,n).*12,...
        0:0.03:3,'LineStyle','none');
    geoshow(land, 'FaceColor',lndcol);
    caxis([0 3]);
    colormap(jet);
    
end

%% Colorbar
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',20);
c=colorbar;
caxis([0 3]);
colormap(jet);
c.Ticks = [0 0.5 1 1.5 2 2.5 3];
c.TickLabels = {'0.0' '0.5' '1.0' '1.5' '2.0' '2.5' '3.0+'};
c.Label.String = 'Interannual variability in {\itF}_{CO2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 20;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_interannual.jpg');
