%% Compare to Hales region

% Index based on Hales et al. (2012) geography
index_Hales = SOCATv2020_grid.lat <= 50 & SOCATv2020_grid.lat >= 22 & ...
        SOCATv2020_grid.distance_from_shore(:,:,1) <= 370;
index_Hales = repmat(index_Hales,1,1,size(SOCATv2020_grid.month_since_1998,1));

% Define seawater pCO2 according to Hales et al. (2012) geography
SOCATv2020_grid.pco2_RF_Hales = SOCATv2020_grid.pco2_RF;
SOCATv2020_grid.pco2_RF_Hales(~index_Hales) = NaN;

% Define atmospheric pCO2 according to Hales et al. (2012) geography
SOCATv2020_grid.pCO2_atm_Hales = SOCATv2020_grid.pCO2_atm;
SOCATv2020_grid.pCO2_atm_Hales(~index_Hales) = NaN;

% Define delta pCO2 according to Hales et al. (2012) geography
SOCATv2020_grid.delta_pCO2_Hales = SOCATv2020_grid.pco2_RF_Hales - ...
                                   SOCATv2020_grid.pCO2_atm_Hales;

% Define CO2 flux according to Hales et al. (2012) geography
SOCATv2020_grid.Fco2_RF_ERA5_Hales = SOCATv2020_grid.Fco2_RF_ERA5;
SOCATv2020_grid.Fco2_RF_ERA5_Hales(~index_Hales) = NaN;

% Define mean pCO2 within Hales et al. (2012) geographic region
SOCATv2020_grid.pco2_RF_Hales_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.pco2_RF_Hales,1,'omitnan'),2,'omitnan'));

% Define mean CO2 flux within Hales et al. (2012) geographic region
SOCATv2020_grid.Fco2_RF_ERA5_Hales_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.Fco2_RF_ERA5_Hales,1,'omitnan'),2,'omitnan'));

% Detrend CO2 flux within Hales et al. (2012) geographic region
[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_Hales_mean,0,2,[6 12]);

% Plot mean CO2 flux over Hales region with trend
figure; hold on;
plot(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_Hales_mean);
plot([SOCATv2020_grid.month_since_1998(1) SOCATv2020_grid.month_since_1998(end)],...
    [x(2).*SOCATv2020_grid.month_since_1998(1)+x(1) ...
    x(2).*SOCATv2020_grid.month_since_1998(end)+x(1)]);
% average over 2001, midpoint of Hales et al. (2012) estimate
mean(x(2).*[37:48]+x(1))

% Plot detrended and deseasonalized mean CO2 flux over Hales region
figure
plot(SOCATv2020_grid.month_since_1998,yr);

% Plot CO2 flux over Hales region
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.delta_pCO2_Hales,3,'omitnan'),...
    -50:1:50,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; colormap(cmocean('thermal',16)); caxis([-50 50]);
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_Hales.jpg');



%% Compare to Turi region

% Index based on Turi et al. (2014) geography
index_Turi = SOCATv2020_grid.lat <= 46 & SOCATv2020_grid.lat >= 33 & ...
        SOCATv2020_grid.distance_from_shore(:,:,1) <= 800;
index_Turi = repmat(index_Turi,1,1,size(SOCATv2020_grid.month_since_1998,1));

% Define seawater pCO2 according to Turi et al. (2014) geography
SOCATv2020_grid.pco2_RF_Turi = SOCATv2020_grid.pco2_RF;
SOCATv2020_grid.pco2_RF_Turi(~index_Turi) = NaN;

% Define atmospheric pCO2 according to Turi et al. (2014) geography
SOCATv2020_grid.pCO2_atm_Turi = SOCATv2020_grid.pCO2_atm;
SOCATv2020_grid.pCO2_atm_Turi(~index_Turi) = NaN;

% Define delta pCO2 according to Turi et al. (2014) geography
SOCATv2020_grid.delta_pCO2_Turi = SOCATv2020_grid.pco2_RF_Turi - ...
                                   SOCATv2020_grid.pCO2_atm_Turi;

% Define CO2 flux according to Turi et al. (2014) geography
SOCATv2020_grid.Fco2_RF_ERA5_Turi = SOCATv2020_grid.Fco2_RF_ERA5;
SOCATv2020_grid.Fco2_RF_ERA5_Turi(~index_Turi) = NaN;

% Define mean pCO2 within Turi et al. (2014) geographic region
SOCATv2020_grid.pco2_RF_Turi_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.pco2_RF_Turi,1,'omitnan'),2,'omitnan'));

% Define mean CO2 flux within Turi et al. (2014) geographic region
SOCATv2020_grid.Fco2_RF_ERA5_Turi_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.Fco2_RF_ERA5_Turi,1,'omitnan'),2,'omitnan'));

% Detrend CO2 flux within Turi et al. (2014) geographic region
[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_Turi_mean,0,2,[6 12]);

% Plot mean CO2 flux over Turi region with trend
figure; hold on;
plot(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_Turi_mean);
plot([SOCATv2020_grid.month_since_1998(1) SOCATv2020_grid.month_since_1998(end)],...
    [x(2).*SOCATv2020_grid.month_since_1998(1)+x(1) ...
    x(2).*SOCATv2020_grid.month_since_1998(end)+x(1)]);
% average over 2001, midpoint of Turi et al. (2014) estimate
mean(x(2).*[25:48]+x(1))

% Plot detrended and deseasonalized mean CO2 flux over Turi region
figure
plot(SOCATv2020_grid.month_since_1998,yr);

% Plot CO2 flux over Turi region
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.latitude(:,:,1),SOCATv2020_grid.longitude(:,:,1),...
    mean(SOCATv2020_grid.Fco2_RF_ERA5_Turi,3,'omitnan'),...
    -3:0.1:3,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([-3 3]); colormap(cmocean('balance','pivot',0));
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
exportgraphics(gcf,'/Users/sharp/Desktop/FCO2_Turi.jpg');



%% Compare to MARCAT CCS region

% Import MARCATs defined by Laruelle et al. (2013)
coscats=shaperead('/Users/sharp/Downloads/hess-17-2029-2013-supplement/GIS Files/Continental_Shelf.shp');

% Index to grid cells that fall within Cal Current MARCAT
index_MARCAT = ...
    inpolygon(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360,coscats(48).Y,coscats(48).X) | ...
    inpolygon(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360,coscats(49).Y,coscats(49).X) | ...
    inpolygon(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360,coscats(50).Y,coscats(50).X) | ...
    inpolygon(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360,coscats(51).Y,coscats(51).X) | ...
    inpolygon(SOCATv2020_grid.lat,SOCATv2020_grid.lon-360,coscats(52).Y,coscats(52).X);
index_MARCAT = repmat(index_MARCAT,1,1,size(SOCATv2020_grid.month_since_1998,1));

% Determine area of each grid cell
SOCATv2020_grid.area = ...
(((SOCATv2020_grid.lat + 0.125) - ...
    (SOCATv2020_grid.lat - 0.125)) .* 110.574) .* ... % latitude distance
(((SOCATv2020_grid.lon + 0.125) - ...
    (SOCATv2020_grid.lon - 0.125)) .* ...
    111.320.*cosd(SOCATv2020_grid.lat)); % longitude distance
% Sum over MARCAT index
SOCATv2020_grid.area_MARCAT = SOCATv2020_grid.area;
SOCATv2020_grid.area_MARCAT(~index_MARCAT(:,:,1)) = NaN;
area_sum_MARCAT = sum(sum(SOCATv2020_grid.area_MARCAT,'omitnan'),'omitnan')./1e3 % Area in 10^3 * km^2

% Define seawater pCO2 according to MARCAT geography
SOCATv2020_grid.pco2_RF_MARCAT = SOCATv2020_grid.pco2_RF;
SOCATv2020_grid.pco2_RF_MARCAT(~index_MARCAT) = NaN;

% Define atmospheric pCO2 according to MARCAT geography
SOCATv2020_grid.pCO2_atm_MARCAT = SOCATv2020_grid.pCO2_atm;
SOCATv2020_grid.pCO2_atm_MARCAT(~index_MARCAT) = NaN;

% Define delta pCO2 according to MARCAT geography
SOCATv2020_grid.delta_pCO2_MARCAT = SOCATv2020_grid.pco2_RF_MARCAT - ...
                                   SOCATv2020_grid.pCO2_atm_MARCAT;

% Define CO2 flux according to MARCAT geography
SOCATv2020_grid.Fco2_RF_ERA5_MARCAT = SOCATv2020_grid.Fco2_RF_ERA5;
SOCATv2020_grid.Fco2_RF_ERA5_MARCAT(~index_MARCAT) = NaN;

% Define mean pCO2 within MARCAT geographic region
SOCATv2020_grid.pco2_RF_MARCAT_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.pco2_RF_MARCAT,1,'omitnan'),2,'omitnan'));

% Define mean CO2 flux within MARCAT geographic region
SOCATv2020_grid.Fco2_RF_ERA5_MARCAT_mean = ...
    squeeze(mean(mean(SOCATv2020_grid.Fco2_RF_ERA5_MARCAT,1,'omitnan'),2,'omitnan'));

% Detrend CO2 flux within MARCAT geographic region
[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_MARCAT_mean,0,2,[6 12]);

% Plot mean over MARCAT region with trend
figure; hold on;
plot(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5_MARCAT_mean);
plot([SOCATv2020_grid.month_since_1998(1) SOCATv2020_grid.month_since_1998(end)],...
    [x(2).*SOCATv2020_grid.month_since_1998(1)+x(1) ...
    x(2).*SOCATv2020_grid.month_since_1998(end)+x(1)]);
% average over 1998-2015, midpoint of MARCAT et al. (2012) estimate
mean(x(2).*[1:18*12]+x(1))
mean(SOCATv2020_grid.Fco2_RF_ERA5_MARCAT_mean(1:18*12),'omitnan')

% Plot detrended and deseasonalized mean over MARCAT region
figure
plot(SOCATv2020_grid.month_since_1998,yr);

% Plot annual mean RF-predicted pCO2
figure; worldmap(latlims,lonlims);
set(gcf,'Position',pos)
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
contourfm(SOCATv2020_grid.lat,SOCATv2020_grid.lon,...
    mean(SOCATv2020_grid.Fco2_RF_ERA5_MARCAT,3,'omitnan'),...
    -3:0.1:3,'LineStyle','none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar; caxis([-3 3]); colormap(cmocean('balance','pivot',0));
c.Label.String = 'pCO_{2} (\muatm)';
c.Label.FontSize = fntsz;
%exportgraphics(gcf,'/Users/sharp/Desktop/pCO2_RF.jpg');




figure;
pcolor(SOCATv2020_grid.lon,SOCATv2020_grid.lat,...
    mean(SOCATv2020_grid.Fco2_RF_ERA5_MARCAT,3,'omitnan'));
colorbar;

figure;
pcolor(SOCATv2020_grid.lon,SOCATv2020_grid.lat,double(index_MARCAT(:,:,1)));

figure; hold on;
mapshow(coscats(48))
mapshow(coscats(49))
mapshow(coscats(50))
mapshow(coscats(51))
mapshow(coscats(52))





