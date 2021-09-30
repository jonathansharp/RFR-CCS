%% This script produces Figure 8 from Sharp et al. (in prep)

% Initialize figure
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol = [1 1 1];
lndcol = [0.5 0.5 0.5];
titlesz = 18;
labelsz = 16;

coast_border = double(SOCATv2021_grid.coast);
coast_border(SOCATv2021_grid.percent_sea < 0.6) = 1;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
tl = axes('Position',[0.1 0.6 0.4 0.35],'Box','on');
tr = axes('Position',[0.56 0.6 0.4 0.35],'Box','on');
bl = axes('Position',[0.05 0.075 0.38 0.45],'Box','on');
br = axes('Position',[0.5 0.075 0.38 0.45],'Box','on');

axes(tl); hold on;
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_coast./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_coast./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
plot(1:12,Fay.Fco2_dom_mean_clim_coast./12,'linewidth',6,'color',[0 0 0],'linestyle',':');
plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_coast./1000,'linewidth',6,'color',[0 0 0.9]);
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('{\itF}CO_{2} (mol C m^{-2} yr^{-1})','fontsize',22);
text(1,-0.18,'coastal','fontsize',18);

axes(tr); hold on;
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_open./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_open./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
p2=plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_open./1000,'linewidth',6,'color',[0 0 0.9]);
p1=plot(1:12,Fay.Fco2_dom_mean_clim_open./12,'linewidth',6,'color',[0 0 0],'linestyle',':');
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
text(1,-0.18,'open ocean','fontsize',18);
legend([p1 p2],{'SeaFlux' 'RFR-CCS'},'location','southeast');

axes(bl); hold on;
worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(Fay.latitude,Fay.longitude,...
    mean(Fay.Fco2_mean,3,'omitnan'),...
    -3:0.5:3,'LineStyle','none');
contourm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    coast_border,'color','k','linewidth',2);
geoshow(land, 'FaceColor',lndcol);
caxis([-3 3]);
colormap(cmocean('balance',20,'pivot',0));
textm(55,233,'SeaFlux','fontsize',26,'color',[0 0 0],'fontweight','bold');

axes(br); hold on;
worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.Fco2_RF_ERA5./1000.*12,3,'omitnan'),...
    -3:0.5:3,'LineStyle','none');
contourm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    coast_border,'color','k','linewidth',2);
geoshow(land, 'FaceColor',lndcol);
caxis([-3 3]);
colormap(cmocean('balance',20,'pivot',0));
textm(55,232,'RFR-CCS','fontsize',26,'color',[0 0 0.9],'fontweight','bold');

% colorbar
cb=axes('Position',[0.05 0.05 0.88 0.5],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
caxis([-3 3]);
colormap(cmocean('balance',12,'pivot',0));
c.Ticks = [-3 0 3];
c.TickLabels = {'-3' '0' '3'};
c.Label.String = '{\itF}CO_{2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 22;
cb.Visible = 'off';

annotation('textbox',[0.11, 0.95, 0, 0],'string','a','fontsize',24)
annotation('textbox',[0.57, 0.95, 0, 0],'string','b','fontsize',24)
annotation('textbox',[0.03, 0.55, 0, 0],'string','c','fontsize',24)
annotation('textbox',[0.48, 0.55, 0, 0],'string','d','fontsize',24)

exportgraphics(gcf,'Figures/Figure8.png');
