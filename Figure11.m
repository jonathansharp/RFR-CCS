%% This script produces Figure 11 from Sharp et al. (in prep)

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

%% Calculate trends in delta pCO2
[trend,pos_trend,neg_trend,no_sea] = ...
gridded_trend(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.delpco2_RF_ERA5,...
    SOCATv2020_grid.pco2_RF,SOCATv2020_grid.lat,y1,y2);

%% Plot trends in delta pCO2
axes(l); worldmap(latlims,lonlims); hold on;
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
s=pcolorm(SOCATv2020_grid.lat-0.125,SOCATv2020_grid.lon-0.125,trend.*12);
set(s, 'EdgeColor', 'none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-5 5]);
colormap(cmocean('balance','pivot',0));
scatterm(SOCATv2020_grid.lat(~pos_trend & ~neg_trend & ~no_sea),...
        SOCATv2020_grid.lon(~pos_trend & ~neg_trend & ~no_sea),...
        '.k');
c.Label.String = 'Trend in \Delta{\itp}CO_{2} (\muatm yr^{-1})';
c.Label.FontSize = fntsz;


%% Calculate trends in delta CO2 flux
[trend,pos_trend,neg_trend,no_sea] = ...
gridded_trend(SOCATv2020_grid.month_since_1998,SOCATv2020_grid.Fco2_RF_ERA5.*12,... % mol C m^-2 yr^-1 month^-1
    SOCATv2020_grid.pco2_RF,SOCATv2020_grid.lat,y1,y2);

%% Plot trends in CO2 flux
axes(r); worldmap(latlims,lonlims); hold on;
setm(gca,'ffacecolor',ocncol);
setm(gca,'fontsize',fntsz);
set(gca,'fontsize',fntsz);
land = shaperead('landareas', 'UseGeoCoords', true);
s=pcolorm(SOCATv2020_grid.lat-0.125,SOCATv2020_grid.lon-0.125,trend.*12.*12.011);
set(s, 'EdgeColor', 'none');
geoshow(land, 'FaceColor',lndcol);
c=colorbar;
caxis([-2 2]);
colormap(cmocean('balance','pivot',0));
scatterm(SOCATv2020_grid.lat(~pos_trend & ~neg_trend & ~no_sea),...
        SOCATv2020_grid.lon(~pos_trend & ~neg_trend & ~no_sea),...
        '.k');
c.Label.String = 'Trend in {\itF}_{CO2} (g C m^{-2} yr^{-2})';
c.Label.FontSize = fntsz;

%% Add figure labels
annotation('textbox',[.03 .85 .1 .1],'String','a','EdgeColor','none','fontsize',32)
annotation('textbox',[.53 .85 .1 .1],'String','b','EdgeColor','none','fontsize',32)
    
%% Export figure
exportgraphics(gcf,'/Users/sharp/Desktop/Figure11.jpg');


%% Embedded function

function [trend,pos_trend,neg_trend,land] = gridded_trend(X,Y,F,L,y1,y2)

% Pre-allocate variables
intercept = nan(size(L));
trend = nan(size(L));
err_intercept = nan(size(L));
err_trend = nan(size(L));
degrees_of_freedom = nan(size(L));
uncert = nan(size(L));

% Determine flux trends for each grid cell
for a = 1:size(L,1)
    for b = 1:size(L,2)
        if sum(~isnan(squeeze(F(a,b,((y1-1998))*12+1:((y2-1998)+1)*12)))) == ((y2-y1)+1)*12
            [yf,yr,x,err,corrmat,r2,n1] = ...
                leastsq2(X(((y1-1998))*12+1:((y2-1998)+1)*12),...
                squeeze(Y(a,b,((y1-1998))*12+1:((y2-1998)+1)*12)),... % 
                0,2,[6 12]);
            intercept(a,b) = x(1);
            trend(a,b) = x(2);
            err_intercept(a,b) = err(1);
            err_trend(a,b) = err(2);
            % Compute autocovariance of residuals
            [acov,acor,lag,dof] = ...
                autocov(X(((y1-1998))*12+1:((y2-1998)+1)*12),yr,36);
            degrees_of_freedom(a,b) = dof;
            %figure; plot(lag,acov);
            % Save scaled uncertainty
            uncert(a,b) = ((err(2)*sqrt(n1)/sqrt(dof-6))*1.68); % 90% confidence
        else
            intercept(a,b) = NaN;
            trend(a,b) = NaN;
            err_intercept(a,b) = NaN;
            err_trend(a,b) = NaN;
            degrees_of_freedom(a,b) = NaN;
            uncert(a,b) = NaN;
        end
    end
end

% Determine if trend is significantly positive or negative
pos_trend = trend > 0 & trend-uncert > 0;
neg_trend = trend < 0 & trend+uncert < 0;
land = isnan(F); land = land(:,:,1);

end