%% 
for i = 1:22
    yrmn(i) = mean(SOCATv2021_grid.pco2_RF_dom_mean((i-1)*12+1:(i-1)*12+12));
end
    
figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 1 0.5]);
set(gca,'fontsize',18);

plot(SOCATv2021_grid.month_since_1998,SOCATv2021_grid.pco2_RF_dom_mean,...
    'Color',rgb('purple'),'linewidth',5);

fill([SOCATv2021_grid.month_since_1998;flipud(SOCATv2021_grid.month_since_1998)],...
    [SOCATv2021_grid.pco2_RF_dom_mean-SOCATv2021_grid.pco2_RF_dom_std';...
    flipud(SOCATv2021_grid.pco2_RF_dom_mean+SOCATv2021_grid.pco2_RF_dom_std')],...
    rgb('purple'),'FaceAlpha',0.25,'LineStyle','none');

%% Figure properties:
xlim([-5 283]); xticks(2*12:5*12:23*12);
ylim([320 460]);
xticklabels({'2000' '2005' '2010' '2015' '2020'});
ylabel('Monthly domain mean {\itp}CO_{2}','fontsize',26);
exportgraphics(gcf,'/Users/sharp/Desktop/domain_mean_pCO2.jpg');

%% 

figure; hold on;
set(gcf,'units','normalized','outerposition',[0.0 0.5 1 0.5]);
set(gca,'fontsize',18);

plot(SOCATv2021_grid.month_since_1998,SOCATv2021_grid.Fco2_RF_dom_mean,...
    'Color',rgb('green'),'linewidth',5);

fill([SOCATv2021_grid.month_since_1998;flipud(SOCATv2021_grid.month_since_1998)],...
    [SOCATv2021_grid.Fco2_RF_dom_mean-SOCATv2021_grid.Fco2_RF_dom_std;...
    flipud(SOCATv2021_grid.Fco2_RF_dom_mean+SOCATv2021_grid.Fco2_RF_dom_std)],...
    rgb('green'),'FaceAlpha',0.25,'LineStyle','none');

%% Figure properties:
xlim([-5 283]); xticks(2*12:5*12:23*12);
%ylim([320 460]);
xticklabels({'2000' '2005' '2010' '2015' '2020'});
ylabel('Monthly domain mean {\itF}_{CO2}','fontsize',26);
exportgraphics(gcf,'/Users/sharp/Desktop/domain_mean_FCO2.jpg');