%% This script produces Figure 5 from Sharp et al. (in prep)

%% Mooring names:
moornames  = {'WA' 'CCE1' 'CCE2' 'NH10' 'LaPush'};
moornames2 = {'Cape Elizabeth' 'CCE1' 'CCE2' 'NH10' 'Chá bâ'};
moornames3 = {'Cape Elizabeth (2006-17)' 'CCE1 (2008-17)' 'CCE2 (2010-17)' 'NH10 (2014-17)' 'Chá bâ (2010-17)'};

% Initialize figure

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol = [1 1 1];
lndcol = [0.94 0.88 0.8];
titlesz = 18;
labelsz = 16;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
tl = axes('Position',[0.1 0.70 0.375 0.25],'Box','on');
tr = axes('Position',[0.6 0.70 0.375 0.25],'Box','on');
ml = axes('Position',[0.1 0.375 0.375 0.25],'Box','on');
mr = axes('Position',[0.6 0.375 0.375 0.25],'Box','on');
bl = axes('Position',[0.1 0.05 0.375 0.25],'Box','on');

aaxx = [tl tr ml mr bl];

% Start loop
for n=1:numel(moornames)
    
    %% Establish axis
    axes(aaxx(n)); hold on;
    
    %% Find closest point on grid to lat-lon
    latidx = ...
        find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)) == ...
        min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{n}).lat)))));
    lonidx = ...
        find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)) == ...
        min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{n}).lon)))));

    %% Determine monthly means for one annual cycle, across mooring available years
    % Pre-allocate
    SOCATv2021_grid.pco2_RF_clim_moor= ...
        nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
    SOCATv2021_grid.pco2_RF_clim_moor_std = ...
        nan(size(SOCATv2021_grid.lon,1),size(SOCATv2021_grid.lat,2),12);
    % Determine monthly means for one annual cycle
    for m = 1:12
        SOCATv2021_grid.pco2_RF_clim_moor(:,:,m) = ...
            mean(SOCATv2021_grid.pco2_RF(:,:,1:12:end),3,'omitnan');
        SOCATv2021_grid.pco2_RF_clim_moor_std(:,:,m) = ...
            std(SOCATv2021_grid.pco2_RF(:,:,1:12:end),0,3,'omitnan');
    end

    %% Plot pCO2 from different sources:
    title(moornames3{n},'FontSize',titlesz);
    set(gca,'fontsize',labelsz);

    fill([1:12 fliplr(1:12)],[[MOORING.(moornames{n}).pCO2SW_monthly(:,1)-...
        mean(MOORING.(moornames{n}).pCO2SW_monthly(:,1)) + ...
        MOORING.(moornames{n}).pCO2SW_monthly(:,2)]' ...
        [flipud(MOORING.(moornames{n}).pCO2SW_monthly(:,1)-...
        mean(MOORING.(moornames{n}).pCO2SW_monthly(:,1)) - ...
        MOORING.(moornames{n}).pCO2SW_monthly(:,2))]'],...
        'k','FaceAlpha',0.25,'LineStyle','none');

    fill([1:12 fliplr(1:12)],[[squeeze(LAR.pCO2_mon(lonidx,latidx,:))-...
        mean(squeeze(LAR.pCO2_mon(lonidx,latidx,:))) + ...
        squeeze(LAR.pCO2_mon_std(lonidx,latidx,:))]' ...
        [flipud(squeeze(LAR.pCO2_mon(lonidx,latidx,:))-...
        mean(squeeze(LAR.pCO2_mon(lonidx,latidx,:))) - ...
        squeeze(LAR.pCO2_mon_std(lonidx,latidx,:)))]'],...
        [0 0.7 0],'FaceAlpha',0.25,'LineStyle','none');
    
    fill([1:12 fliplr(1:12)],[[squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))-...
        mean(squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))) + ...
        squeeze(SOCATv2021_grid.pco2_RF_clim_std_2015(lonidx,latidx,:))]' ...
        [flipud(squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))-...
        mean(squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:)))) - ...
        squeeze(SOCATv2021_grid.pco2_RF_clim_std_2015(lonidx,latidx,:))]'],...
        [0 0 0.9],'FaceAlpha',0.25,'LineStyle','none');
    
    zer=plot([0 13],[0 0],':','LineWidth',3,'Color','k');
    
    mo=plot(1:12,MOORING.(moornames{n}).pCO2SW_monthly(:,1)-...
        mean(MOORING.(moornames{n}).pCO2SW_monthly(:,1)),'-o',...
        'LineWidth',3,'Color','k');
    mo1=scatter(1:12,MOORING.(moornames{n}).pCO2SW_monthly(:,1)-...
        mean(MOORING.(moornames{n}).pCO2SW_monthly(:,1)),80,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'LineWidth',1);

    land=plot(1:12,squeeze(LAND.pCO2(lonidx,latidx,:))-...
        mean(squeeze(LAND.pCO2(lonidx,latidx,:))),'-o',...
        'LineWidth',3,'Color',[.5 0 .5]);
    scatter(1:12,squeeze(LAND.pCO2(lonidx,latidx,:))-...
        mean(squeeze(LAND.pCO2(lonidx,latidx,:))),80,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[.5 0 .5],...
        'LineWidth',1);
    
    lar=plot(1:12,squeeze(LAR.pCO2_mon(lonidx,latidx,:))-...
        mean(squeeze(LAR.pCO2_mon(lonidx,latidx,:))),'-o',...
        'LineWidth',3,'Color',[0 0.7 0]);
    scatter(1:12,squeeze(LAR.pCO2_mon(lonidx,latidx,:))-...
        mean(squeeze(LAR.pCO2_mon(lonidx,latidx,:))),80,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 0.7 0],...
        'LineWidth',1);
    
    rf=plot(1:12,squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))-...
        mean(squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))),'-o',...
        'LineWidth',3,'Color',[0 0 0.9]);
    scatter(1:12,squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))-...
        mean(squeeze(SOCATv2021_grid.pco2_RF_clim_2015(lonidx,latidx,:))),80,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.9],...
        'LineWidth',1);
    
    %% Figure properties:
    xlim([0.5 12.5]);
    ylim([-150 150]);
    ylabel('pCO_{2} Anomaly (\muatm)');
    xticks([1:12]);
    xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
end

%% Legends for figures:
br = axes('Position',[0.6 0.05 0.375 0.25],'Box','off');
hl=legend(br,[mo land lar rf],{'Mooring' 'L20' 'L17' 'RFR-CCS-clim'});
set(hl,'Position', [0.65 0.075 0.225 0.2],'fontsize',28);
%l = findobj(l,'Type','Line');
%set(l,'MarkerSize',10);
br.Visible = 'off';

%% Export figure
exportgraphics(gcf,['Figures/Figure5.png']);
