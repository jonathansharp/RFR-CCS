% Monthly CO2 flux standard deviation
%
% This function determines the monthly standard deviation of CO2 flux for
% individual grid cells of RFR-NEP across the timeframe of 1998 - 2019
% This gives indication of the interannual variability in CO2 flux

for m=1:12
    SOCATv2020_grid.Fco2_monthly_std(:,:,m) = ...
        std(SOCATv2020_grid.Fco2_RF_ERA5(:,:,m:12:end),[],3);
end

% Monthly Mean Standard Deviations in FCO2 in Space
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.Fco2_monthly_std(:,:,:),3,'omitnan'),...
    0:0.1:3,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
colormap(cmocean('thermal'));
caxis([0 3]);
%c.TickLabels = {'300' '320' '340' '360' '380' '400' '420' '440+'};
c.Label.String = 'Monthly Mean Standard Deviation in FCO2';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_std_spatial.jpg');

% Spatial Mean Standard Deviations in FCO2 by Month
figure;
set(gca,'fontsize',16);
s1=scatter(1:12,squeeze(mean(mean(SOCATv2020_grid.Fco2_RF_ERA5_avg(:,:,:),1,'omitnan'),2,'omitnan')),...
    50,'ok','filled'); hold on;
e1=errorbar(1:12,squeeze(mean(mean(SOCATv2020_grid.Fco2_RF_ERA5_avg(:,:,:),1,'omitnan'),2,'omitnan')),...
    squeeze(mean(mean(SOCATv2020_grid.Fco2_monthly_std(:,:,:),1,'omitnan'),2,'omitnan')),'-k',...
    'linewidth',2);
xlim([0 13]); xticks([1:12]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
xlabel('Month of Year');
ylabel('Domain monthly mean FCO_{2} (mol m^{-2} yr^{-1}) from 1998-2019');
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_std_temporal.jpg');
