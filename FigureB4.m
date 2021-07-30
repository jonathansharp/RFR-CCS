%% This script produces Figure B4 from Sharp et al. (in prep)

%% Determine correlation between temperature and pCO2
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2021_grid.pco2_RF(a,b,:)))) == 0
            [SOCATv2021_grid.t_vs_pco2_corr(a,b,:),...
             SOCATv2021_grid.p_t_corr(a,b,:)] = ...
                corr(squeeze(SOCATv2021_grid.SST(a,b,1:276)),...
                     squeeze(SOCATv2021_grid.pco2_RF(a,b,1:276)));
        else
            SOCATv2021_grid.t_vs_pco2_corr(a,b,:) = NaN;
            SOCATv2021_grid.p_t_corr(a,b,:) = NaN;
        end
    end
end

%% Determine correlation between winds and pCO2
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2021_grid.pco2_RF(a,b,:)))) == 0
            [SOCATv2021_grid.wind_vs_pco2_corr(a,b,:),...
             SOCATv2021_grid.p_wind_corr(a,b,:)] = ...
                corr(squeeze(SOCATv2021_grid.wind_speed(a,b,1:276)),...
                     squeeze(SOCATv2021_grid.pco2_RF(a,b,1:276)));
        else
            SOCATv2021_grid.wind_vs_pco2_corr(a,b,:) = NaN;
            SOCATv2021_grid.p_wind_corr(a,b,:) = NaN;
        end
    end
end

%% Initialize figure
figure;
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
tl = axes('Position',[0.02 0.55 0.425 0.4],'Box','on');
tr = axes('Position',[0.42 0.55 0.425 0.4],'Box','on');
bl = axes('Position',[0.02 0.05 0.425 0.4],'Box','on');
br = axes('Position',[0.42 0.05 0.425 0.4],'Box','on');

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
ocncol = [1 1 1];
lndcol = [1 1 1];
fntsz = 18;

%% Plot correlation between temperature and pCO2 (Xa)
axes(tl); title('Temp. vs. {\itp}CO_{2(sw)}'); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    SOCATv2021_grid.t_vs_pco2_corr);
geoshow(land,'FaceColor',lndcol,'linewidth',1);
caxis(tl,[-1 1]);
colormap(tl,cmocean('balance','pivot',0));% c.Label.String = 'Correlation (MLD vs. pCO_{2})';
% c.Label.FontSize = 18;
% exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr.jpg');

%% Plot p-value of correlation between temperature and pCO2 (Xb)
axes(bl); title('{\itp}-value'); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    SOCATv2021_grid.p_t_corr);
geoshow(land,'FaceColor',lndcol,'linewidth',1);
caxis(bl,[0 0.5]);
colormap(bl,cmocean('balance','pivot',0.05));% c.Label.String = 'Correlation (MLD vs. pCO_{2})';
% c.Label.FontSize = 18;
% exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr.jpg');

%% Plot correlation between wind speed and pCO2 (Xc)
axes(tr); title('Wind Speed vs. {\itp}CO_{2(sw)}'); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    SOCATv2021_grid.wind_vs_pco2_corr);
geoshow(land,'FaceColor',lndcol,'linewidth',1);
caxis(tr,[-1 1]);
colormap(tr,cmocean('balance','pivot',0));% c.Label.String = 'Correlation (MLD vs. pCO_{2})';
% c.Label.FontSize = 18;
% exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr.jpg');


%% Plot p-value of correlation between wind speed and pCO2 (Xd)
axes(br); title('{\itp}-value'); worldmap(latlims,lonlims);
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',14);
set(gca,'fontsize',14);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,...
    SOCATv2021_grid.p_wind_corr);
geoshow(land,'FaceColor',lndcol,'linewidth',1);
caxis(br,[0 0.5]);
colormap(br,cmocean('balance','pivot',0.05));% c.Label.String = 'Correlation (MLD vs. pCO_{2})';
% c.Label.FontSize = 18;
% exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr.jpg');



%% Colorbar
ax=axes('Position',[0.45 0.52 0.45 0.45],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
caxis(ax,[-1 1]);
colormap(ax,cmocean('balance','pivot',0));
% c.Ticks = [-50 -25 0 25 50];
% c.TickLabels = {'-50' '-25' '0' '25' '50+'};
c.Label.String = 'Correlation Coefficient';
c.Label.FontSize = 24;
ax.Visible = 'off';

ax=axes('Position',[0.45 0.02 0.45 0.45],'Box','off');
set(gca,'fontsize',24);
c=colorbar;
caxis(ax,[0 0.5]);
colormap(ax,cmocean('balance','pivot',0.05));% c.Ticks = [-50 -25 0 25 50];
% c.TickLabels = {'-50' '-25' '0' '25' '50+'};
c.Label.String = '{\itp}-value';
c.Label.FontSize = 24;
ax.Visible = 'off';

%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure4B.jpg');
