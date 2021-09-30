% Plots SOCAT observations by month

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol  = [1 1 1];
lndcol = [0.5 0.5 0.5];
titlesz = 20;
labelsz = 22;
fontsz = 12;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.675 1]);

mnth = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
x_ax = [0.02 0.23 0.43 0.63 0.02 0.23 0.43 0.63 0.02 0.23 0.43 0.63];
y_ax = [0.71 0.71 0.71 0.71 0.38 0.38 0.38 0.38 0.05 0.05 0.05 0.05];

for n = 1:12

    % Define axis
    axes('Position',[x_ax(n) y_ax(n) 0.25 0.25],'Box','off');
    
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

exportgraphics(gcf,'/Users/sharp/Desktop/FigureB1.png');
