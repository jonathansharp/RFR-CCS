%% This script produces Figure 2 from Sharp et al. (in prep)
% It plots pCO2 at the CCE2 mooring, along with pCO2 from corresponding
% grid cells in RFR-CCS, RFR-CCS-Eval, and Laruelle et al. (2017).

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [0 0 1 0.5];
titlesz = 28;
fontsz = 22;

% Match latitude in climatology
latidx = ...
    find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.CCE2.lat)) == ...
    min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.CCE2.lat)))));
lonidx = ...
    find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.CCE2.lon)) == ...
    min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.CCE2.lon)))));

% Determine mooring month since 1998
MOORING.CCE2.month_since_1998 = ...
    (MOORING.CCE2.date(:,1)-1998).*12 + MOORING.CCE2.date(:,2);

% Determine means of each month since 1998
MOORING.CCE2.pCO2SW_mon = nan(max(SOCATv2021_grid.month_since_1998),1);
for m=1:max(SOCATv2021_grid.month_since_1998)
    MOORING.CCE2.month_since_1998_mon(m) = m;
    MOORING.CCE2.pCO2SW_mon_mean(m) = ...
        mean(MOORING.CCE2.pCO2SW(MOORING.CCE2.month_since_1998==m),'omitnan');
    MOORING.CCE2.pCO2SW_mon_std(m)  = ...
        std(MOORING.CCE2.pCO2SW(MOORING.CCE2.month_since_1998==m),'omitnan');
end

% Plot climatology pCO2 and mooring measured pCO2
figure; box on;
set(gcf,'units','normalized','outerposition',pos);
set(gca,'fontsize',fontsz);
title('CCE2','fontsize',titlesz);
hold on
date = datenum([repmat(1998,276,1) SOCATv2021_grid.month_since_1998 ones(276,1)]);
p1=plot(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)),'linewidth',4,'Color',[0 0 0.9]);
s1=scatter(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)),100,'^',...
    'MarkerFaceColor',[0 0 0.9],'MarkerEdgeColor','k');
p2=plot(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)),':','linewidth',4,'Color',[0 0 0.9]);
s2=scatter(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)),100,...
    'MarkerFaceColor',[0 0 0.9],'MarkerEdgeColor','k');
date = datenum([repmat(1998,276,1) MOORING.CCE2.month_since_1998_mon' ones(276,1)]);
p3=plot(date,MOORING.CCE2.pCO2SW_mon_mean,'k','linewidth',4);
s3=scatter(date,MOORING.CCE2.pCO2SW_mon_mean,100,...
    'MarkerFaceColor','k','MarkerEdgeColor','k');
p4=plot(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)),'Color',[0 0.7 0],'linewidth',4);
s4=scatter(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)),100,...
    'MarkerFaceColor',[0 0.7 0],'MarkerEdgeColor','k');
legend([p3 p1 p2 p4],{'Mooring Observations'...
                       'RFR-CCS'...
                       'RFR-CCS-Eval'...
                       'L17'},'location','northeast','box','on');
ylabel('{\itp}CO_{2(sw)} (\muatm)');
xlabel('Year');
datetick('x','yyyy');
xlim([min(date(~isnan(MOORING.CCE2.pCO2SW_mon_mean))) ...
      max(date(~isnan(MOORING.CCE2.pCO2SW_mon_mean)))]);
ylim([300 550]);

%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure2.jpg');

