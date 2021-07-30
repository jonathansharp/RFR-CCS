% Plots SOCAT observations by month

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol  = [0.94 0.97 1.0];
lndcol  = [1 1 1];
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
    title(mnth{n},'fontsize',titlesz);
    setm(gca,'fontsize',fontsz);
    set(gca,'fontsize',fontsz);
    land = shaperead('landareas', 'UseGeoCoords', true);
    pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
        (sum(~isnan(SOCATv2021_grid.all.pco2_ave_weighted(:,:,n:12:end)),3)));
    caxis([0 6]);
    myColorMap = jet(6);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    geoshow(land, 'FaceColor',lndcol,'linewidth',1);

end

% Define axis
ax=axes('Position',[0.8 0.05 0.15 0.9],'Box','off');
set(gca,'fontsize',labelsz);
myColorMap = jet(6);
myColorMap(1,:) = 1;
colormap(myColorMap);
c=colorbar;
caxis([0 6]);
c.Ticks = [0.5 1.5 2.5 3.5 4.5 5.5];
c.TickLabels = {'0' '1' '2' '3' '4' '5+'};
c.Label.String = 'Number of Years Represented';
c.Label.FontSize = labelsz;
ax.Visible = 'off';

exportgraphics(gcf,'/Users/sharp/Desktop/FigureB1.jpg');
