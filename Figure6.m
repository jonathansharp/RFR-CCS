%% This script produces Figure 6 from Sharp et al. (in prep)

%% Initialize figure

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [617, 599, 820, 820];
ocncol  = [1 1 1];
lndcol  = [0.5 0.5 0.5];
titlesz = 28;
labelsz = 22;
fontsz = 18;

figure;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

l = axes('Position',[0.05 0.05 0.4 0.9],'Box','on');
r = axes('Position',[0.55 0.05 0.4 0.9],'Box','on');

%% Plot annual mean RF-predicted pCO2
axes(l); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fontsz);
set(gca,'fontsize',fontsz);
title('Annual Mean','fontsize',titlesz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    mean(SOCATv2021_grid.pco2_RF_clim,3,'omitnan'),...
    300:10:440,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol,'linestyle','none');
moornames={'CCE1' 'CCE2' 'NH10' 'WA' 'LaPush'};
moornames2={'CCE1' 'CCE2' 'NH10' 'Cape Elizabeth' 'Châ bá'};
for n=1:numel(moornames)
    % Find closest point in RF climatology to lat-lon
    latidx = ...
        find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    scatterm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1),...
        400,mean(MOORING.(moornames{n}).pCO2SW,'omitnan'),'filled',...
        'o','MarkerEdgeColor','w','LineWidth',2);
    if n==5 % Cha ba
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==3 % NH10
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+2,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==1 % CCE1
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+5.5,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==2 % CCE2
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+0.5,SOCATv2021_grid.longitude(lonidx,latidx,1)+1.5,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    else
        textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    end
end
c=colorbar;
%cptcmap('GMT_ocean','flip',true);
colormap(l,cmocean('haline',14));
caxis([300 440]);
c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440'};
c.Label.String = '{\itp}CO_{2(sw)} (\muatm)';
c.Label.FontSize = labelsz;

%% Plot RF-predicted pCO2 amplitude
axes(r); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fontsz);
set(gca,'fontsize',fontsz);
title('Seasonal Amplitude','fontsize',titlesz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    log10(SOCATv2021_grid.pco2_RF_amp),...
    1:0.1:2.5,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol,'linestyle','none');
for n=1:numel(moornames)
    % Find closest point in RF climatology to lat-lon
    latidx = ...
        find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));
    scatterm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1),...
        400,log10(MOORING.(moornames{n}).pCO2SW_amplitude),'filled',...
        'MarkerEdgeColor','k','LineWidth',2);
    if n==5 % Cha ba
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+1,SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==3 % NH10
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+2,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==1 % CCE1
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+5.5,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    elseif n==2 % CCE2
    textm(SOCATv2021_grid.latitude(lonidx,latidx,1)+0.5,SOCATv2021_grid.longitude(lonidx,latidx,1)+1.5,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    else
        textm(SOCATv2021_grid.latitude(lonidx,latidx,1),SOCATv2021_grid.longitude(lonidx,latidx,1)+3,...
        moornames2{n},'Color','w','FontWeight','bold','fontsize',labelsz);
    end
end
c=colorbar;
caxis([1 2.5]);
colormap(r,cmocean('thermal',15));
c.Ticks = [1 1.477 2 2.477];
c.TickLabels = [10 30 100 300];
c.Label.String = '{\itp}CO_{2(sw)} (\muatm)';
c.Label.FontSize = labelsz;

%% Add figure labels
annotation('textbox',[.03 .85 .1 .1],'String','a','EdgeColor','none','fontsize',32)
annotation('textbox',[.53 .85 .1 .1],'String','b','EdgeColor','none','fontsize',32)
    
%% Export figure
exportgraphics(gcf,'Figures/Figure6.png');

