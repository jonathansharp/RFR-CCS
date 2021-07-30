%% This script produces Figure 9 from Sharp et al. (in prep)

%% Initialize figure
figure;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
l = axes('Position',[0.05 0.05 0.4 0.9],'Box','on');
r = axes('Position',[0.55 0.05 0.4 0.9],'Box','on');

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol  = [1 1 1];
lndcol  = [1 1 1];
fntsz   = 18;

%% Plot annual mean RF-predicted pCO2
axes(l); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
title('Annual Mean','fontsize',28);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    mean(SOCATv2021_grid.Fco2_RF_ERA5.*12,3,'omitnan'),...
    -3:0.3:3,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
% moornames={'CCE1' 'CCE2' 'NH10' 'WA' 'LaPush'};
% moornames2={'CCE1' 'CCE2' 'NH10' 'Cape Elizabeth' 'Châ bá'};
% for n=1:numel(moornames)
%     % Find closest point in RF climatology to lat-lon
%     latidx = ...
%         find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
%         min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
%     lonidx = ...
%         find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
%         min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
%     scatterm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1),...
%         400,mean(MOORING.(moornames{n}).pCO2SW,'omitnan'),'filled',...
%         'o','MarkerEdgeColor','k','LineWidth',2);
%     if n==6
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
%     elseif n==1
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+6,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
%     elseif n==2
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+2,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
%     else
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
%     end
% end
c=colorbar;
caxis([-3 3]);
colormap(cmocean('balance',20,'pivot',0));
c.Ticks = [-3 -2 -1 0 1 2 3];
c.TickLabels = {'-3' '-2' '-1' '0' '1' '2' '3'};
c.Label.String = '{\itF}_{CO2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 20;

%% Plot RF-predicted pCO2 amplitude
axes(r); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
title('Seasonal Amplitude','fontsize',22);
contourfm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    log10(SOCATv2021_grid.Fco2_RF_ERA5_amp),...
    -1.3:0.05:-0.1,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
% for n=1:numel(moornames)
%     % Find closest point in RF climatology to lat-lon
%     latidx = ...
%         find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
%         min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
%     lonidx = ...
%         find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
%         min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
%     scatterm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1),...
%         400,log10(MOORING.(moornames{n}).pCO2SW_amplitude),'filled',...
%         'MarkerEdgeColor','k','LineWidth',2);
%     if n==6
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
%     elseif n==1
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+6,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
%     elseif n==2
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+2,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16); 
%     else
%     textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
%         moornames2{n},'Color','k','FontWeight','bold','fontsize',16);
%     end
% end
c=colorbar;
caxis([-1.3 -0.1]);
colormap(r,cmocean('thermal',24));
c.Ticks = [-1.3 -1 -0.69897 -0.39794 -0.1549];
c.TickLabels = [0.05 0.1 0.2 0.4 0.7];
c.Label.String = '{\itF}_{CO2} (mol C m^{-2} yr^{-1})';
c.Label.FontSize = 20;

%% Add figure labels
annotation('textbox',[.03 .85 .1 .1],'String','a','EdgeColor','none','fontsize',32)
annotation('textbox',[.53 .85 .1 .1],'String','b','EdgeColor','none','fontsize',32)
    
%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure9.jpg');

