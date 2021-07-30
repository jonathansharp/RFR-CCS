%% This script produces Figure A3 from Sharp et al. (in prep)

latlims = [latmin latmax];
lonlims = [lonmin lonmax];
pos = [0 0 1 0.5];
titlesz = 28;
fontsz = 22;


% Importing chlorophyll
path = '/Volumes/2TB Hard Drive/SATELLITE_DATA/';
load(strcat(path,'CHL.mat')); clear path
CHL.latitude = repmat(CHL.lat',1,size(CHL.chl,2),size(CHL.chl,3));
CHL.longitude = repmat(CHL.lon,size(CHL.chl,1),1,size(CHL.chl,3));
% Cut down dataset to CCS limits
lonidx = CHL.lon >= lonmin & CHL.lon <= lonmax;
latidx = CHL.lat >= latmin & CHL.lat <= latmax;
CHL.chl = CHL.chl(latidx,lonidx,:);
CHL.latitude = CHL.latitude(latidx,lonidx,:);
CHL.longitude = CHL.longitude(latidx,lonidx,:);
CHL.longitude = CHL.longitude + 360; % Convert to 360 degrees longitude
% Match time frame of SOCAT data
CHL.date = datevec(CHL.time);
CHL.month_since_1998 = (CHL.date(:,1)-1998).*12 + CHL.date(:,2);
idx = CHL.month_since_1998 >= min(SOCATv2021_grid.month_since_1998) & ...
      CHL.month_since_1998 <= max(SOCATv2021_grid.month_since_1998);
CHL.chl = CHL.chl(:,:,idx);
CHL.latitude = CHL.latitude(:,:,idx);
CHL.longitude = CHL.longitude(:,:,idx);
% Interpolate onto SOCAT grid
for t = 1:max(SOCATv2021_grid.month_since_1998)
    interp = griddedInterpolant(flipud(CHL.longitude(:,:,t))',flipud(CHL.latitude(:,:,t))',flipud(CHL.chl(:,:,t))');
    supp_fig.CHL(:,:,t) = interp(SOCATv2021_grid.longitude(:,:,1),SOCATv2021_grid.latitude(:,:,1));
end
% Eliminate lake values
supp_fig.CHL(SOCATv2021_grid.latitude > 40 & ...
    SOCATv2021_grid.longitude(:,:,1) > 240) = NaN;
% Interpolate over some gaps in CHL dataset (linear, 1-D, time)
for g = 1:size(supp_fig.CHL,1)
    for h = 1:size(supp_fig.CHL,2)
        if sum(~isnan(supp_fig.CHL(g,h,:))) >= 100
            Chl = squeeze(supp_fig.CHL(g,h,:));
            idx = ~isnan(Chl);
            Chlfit = interp1(SOCATv2021_grid.month_since_1998(idx),...
                Chl(idx),SOCATv2021_grid.month_since_1998,'linear');
            supp_fig.CHL2(g,h,:) = Chlfit;
        else
            supp_fig.CHL2(g,h,:) = NaN;
        end
    end
end
% Interpolate over remaining gaps in CHL dataset (nearest, 1-D, time)
for g = 1:size(supp_fig.CHL2,1)
    for h = 1:size(supp_fig.CHL2,2)
        if sum(~isnan(supp_fig.CHL2(g,h,:))) >= 100
            Chl = squeeze(supp_fig.CHL2(g,h,:));
            idx = ~isnan(Chl);
            Chlfit = interp1(SOCATv2021_grid.month_since_1998(idx),...
                Chl(idx),SOCATv2021_grid.month_since_1998,'nearest','extrap');
            supp_fig.CHL2(g,h,:) = Chlfit;
        else
            supp_fig.CHL2(g,h,:) = NaN;
        end
    end
end
% Replace values less than zero with value close to zero
supp_fig.CHL2(supp_fig.CHL2<0) = 0.001;

%
figure; box on;
set(gcf,'units','normalized','outerposition',pos);
set(gca,'fontsize',fontsz);
hold on
p1=plot(SOCATv2021_grid.month_since_1998,squeeze(supp_fig.CHL2(20,160,:)),'s-','linewidth',3)
p2=plot(SOCATv2021_grid.month_since_1998,squeeze(supp_fig.CHL(20,160,:)),'s-','linewidth',3)
legend([p2 p1],{'CHL observations' 'Interpolated points'},'location','northwest','box','on');
ylabel('CHL (mg m^{-2})','fontsize',fontsz);
xlim([min(SOCATv2021_grid.month_since_1998) ...
      max(SOCATv2021_grid.month_since_1998)]);
xticks(2*12:5*12:22*12);
xticklabels({'2000' '2005' '2010' '2015' '2020'});

% Save figure
exportgraphics(gcf,'/Users/sharp/Desktop/FigureA3.jpg');

%clear supp_fig