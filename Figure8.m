%% This script produces Figure 8 from Sharp et al. (in prep)

% Initialize figure
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol = [1 1 1];
lndcol = [0.5 0.5 0.5];
titlesz = 20;
labelsz = 18;

coast_border = double(SOCATv2021_grid.coast);
coast_border(SOCATv2021_grid.percent_sea < 0.6) = 1;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.9 1])
tl = axes('Position',[0.1 0.6 0.4 0.35],'Box','on');
tr = axes('Position',[0.53 0.6 0.4 0.35],'Box','on');
bl = axes('Position',[0.05 0.075 0.28 0.45],'Box','on');
bm = axes('Position',[0.3 0.075 0.28 0.45],'Box','on');
br = axes('Position',[0.65 0.075 0.28 0.45],'Box','on');

std_coast = std([Fay.Fco2_CMEMS_FFNN_dom_mean_clim_coast;...
                 Fay.Fco2_CSIR_ML6_dom_mean_clim_coast;...
                 Fay.Fco2_JENA_MLS_dom_mean_clim_coast;...
                 Fay.Fco2_JMA_MLR_dom_mean_clim_coast;...
                 Fay.Fco2_MPI_SOMFFN_dom_mean_clim_coast;...
                 Fay.Fco2_NIES_FNN_dom_mean_clim_coast]);
std_open = std([Fay.Fco2_CMEMS_FFNN_dom_mean_clim_open;...
                 Fay.Fco2_CSIR_ML6_dom_mean_clim_open;...
                 Fay.Fco2_JENA_MLS_dom_mean_clim_open;...
                 Fay.Fco2_JMA_MLR_dom_mean_clim_open;...
                 Fay.Fco2_MPI_SOMFFN_dom_mean_clim_open;...
                 Fay.Fco2_NIES_FNN_dom_mean_clim_open]);

axes(tl); hold on;
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_coast./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_coast./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
% fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_clim_coast./12 + std_coast./12....
%     fliplr(Fay.Fco2_dom_mean_minus_clim_coast./12 - std_coast./12)],...
%     [0 0 0],'linestyle','none','FaceAlpha',0.2);
p1=plot(1:12,Fay.Fco2_dom_mean_clim_coast./12,'linewidth',6,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_CMEMS_FFNN_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_CSIR_ML6_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_JENA_MLS_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_JMA_MLR_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_MPI_SOMFFN_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_NIES_FNN_dom_mean_clim_coast./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
p2=plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_coast./1000,'linewidth',6,'color',clr1);
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
ylabel('{\itF}CO_{2} (mol C m^{-2} yr^{-1})','fontsize',22);
text(1,-0.18,'coastal','fontsize',22);
legend([p1 p2],{'SeaFlux' 'RFR-CCS'},'location','southeast');

axes(tr); hold on;
set(gca,'fontsize',18);
fill([1:12 12:-1:1],[Fay.Fco2_dom_mean_plus_clim_open./12 ...
    fliplr(Fay.Fco2_dom_mean_minus_clim_open./12)],...
    [0 0 0],'linestyle','none','FaceAlpha',0.2);
p2=plot(1:12,SOCATv2021_grid.Fco2_RF_dom_mean_clim_open./1000,'linewidth',6,'color',clr1);
plot(1:12,Fay.Fco2_CMEMS_FFNN_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_CSIR_ML6_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_JENA_MLS_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_JMA_MLR_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_MPI_SOMFFN_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
plot(1:12,Fay.Fco2_NIES_FNN_dom_mean_clim_open./12,'linewidth',1,'color',[0 0 0],'linestyle',':');
p1=plot(1:12,Fay.Fco2_dom_mean_clim_open./12,'linewidth',6,'color',[0 0 0],'linestyle',':');
xlim([0 13]); xticks([1:12]); ylim([-0.2 0.1]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
text(1,-0.18,'open ocean','fontsize',22);
legend([p1 p2],{'SeaFlux' 'RFR-CCS'},'location','southeast');

axes(bl); hold on;
worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.Fco2_RF_ERA5./1000.*12,3,'omitnan'),...
    -3.25:0.5:3.25,'LineStyle','none');
contourm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    coast_border,'color','k','linewidth',2);
geoshow(land, 'FaceColor',lndcol);
set(bl,'CLim',[-3.25 3.25]);
set(bl,'colormap',cmocean('balance',13,'pivot',0));
textm(55,232,'RFR-CCS','fontsize',26,'color',clr1,'fontweight','bold');

axes(bm); hold on;
worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(Fay.latitude,Fay.longitude,...
    mean(Fay.Fco2_mean,3,'omitnan'),...
    -3.25:0.5:3.25,'LineStyle','none');
contourm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    coast_border,'color','k','linewidth',2);
geoshow(land, 'FaceColor',lndcol);
set(bm,'CLim',[-3.25 3.25]);
set(bm,'colormap',cmocean('balance',13,'pivot',0));
textm(55,233,'SeaFlux','fontsize',26,'color',[0 0 0],'fontweight','bold');

% colorbar
cb1=axes('Position',[0.05 0.05 0.55 0.5],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
caxis([-3 3]);
set(cb1,'CLim',[-3.25 3.25]);
set(cb1,'colormap',cmocean('balance',13,'pivot',0));
c.Ticks = [-3 0 3];
c.TickLabels = {'-3' '0' '3'};
c.Label.String = '{\itF}CO_{2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 22;
cb1.Visible = 'off';

axes(br); hold on;
worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
temp_mean = interp2(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    mean(SOCATv2021_grid.Fco2_RF_ERA5./1000.*12,3,'omitnan'),...
    Fay.latitude,Fay.longitude+360);
contourfm(Fay.latitude,Fay.longitude,temp_mean-...
    mean(Fay.Fco2_mean,3,'omitnan'),-1.0625:0.125:1.0625,'LineStyle','none');
contourm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    coast_border,'color','k','linewidth',2);
geoshow(land, 'FaceColor',lndcol);
set(br,'CLim',[-1.0625 1.0625]);
set(br,'colormap',cmocean('balance',17,'pivot',0));
textm(55,233,'RFR-CCS','fontsize',20,'color',clr1,'fontweight','bold');
textm(53,237,'– SeaFlux','fontsize',20,'color',[0 0 0],'fontweight','bold');

% colorbar
cb2=axes('Position',[0.05 0.05 0.9 0.5],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
set(cb2,'CLim',[-1.0625 1.0625]);
set(cb2,'colormap',cmocean('balance',17,'pivot',0));
c.Ticks = [-1 0 1];
c.TickLabels = {'-1' '0' '1'};
c.Label.String = '\Delta{\itF}CO_{2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 22;
cb2.Visible = 'off';

annotation('textbox',[0.11, 0.95, 0, 0],'string','a','fontsize',26)
annotation('textbox',[0.54, 0.95, 0, 0],'string','b','fontsize',26)
annotation('textbox',[0.05, 0.55, 0, 0],'string','c','fontsize',26)
annotation('textbox',[0.30, 0.55, 0, 0],'string','d','fontsize',26)
annotation('textbox',[0.66, 0.55, 0, 0],'string','e','fontsize',26)

exportgraphics(gcf,'Figures/Figure8.png');
