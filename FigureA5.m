%% This script produces Figure 7 from Sharp et al. (in prep)
% It gives monthly mean fields of pCO2 determine by a random forest
% regression on SOCATv2021 observations in the northeast Pacific

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol = [1 1 1];
lndcol = [0.5 0.5 0.5];
titlesz = 20;
labelsz = 22;
fontsz = 12;

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
    setm(gca,'ffacecolor',ocncol);
    setm(gca,'fontsize',fontsz);
    set(gca,'fontsize',fontsz);
    title(mnth{n},'fontsize',titlesz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    contourfm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
        mean(SOCATv2021_grid.SST(:,:,n:12:end),3,'omitnan'),...
        0:5:35,'LineStyle','none');
    geoshow(land, 'FaceColor',lndcol,'linewidth',1);
    caxis([0 35]);
    colormap(cmocean('haline',14));
    
end

%% Colorbar
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',labelsz);
c=colorbar;
caxis([0 35]);
colormap(cmocean('haline',14));
%c.Ticks = [300 320 340 360 380 400 420 440];
%c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'SST (degC)';
c.Label.FontSize = labelsz;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/FigureA5.jpg');
