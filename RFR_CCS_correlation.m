% Determine correlation between temperature and pCO2
for a=1:140
    for b = 1:180
        if sum(isnan(squeeze(SOCATv2021_grid.pco2_RF(a,b,:)))) == 0
            [SOCATv2021_grid.t_vs_pco2_corr(a,b,:),...
             SOCATv2021_grid.p_corr(a,b,:)] = ...
                corr(squeeze(SOCATv2021_grid.MLD(a,b,1:276)),...
                     squeeze(SOCATv2021_grid.pco2_RF(a,b,1:276)));
        else
            SOCATv2021_grid.t_vs_pco2_corr(a,b,:) = NaN;
            SOCATv2021_grid.p_corr(a,b,:) = NaN;
        end
    end
end

% Plot temperature vs. pCO2 correlation
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
figure; worldmap(latlims,lonlims);
set(gcf,'Position',[617, 599, 820, 820])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',18);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    SOCATv2021_grid.t_vs_pco2_corr);
geoshow(land, 'FaceColor','w');
c=colorbar; caxis([-1 1]); colormap(cmocean('balance','pivot',0));
c.Label.String = 'Correlation (MLD vs. pCO_{2})';
c.Label.FontSize = 18;
exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr.jpg');

% Plot p-values for temperature vs. pCO2 correlation
latlims = [latmin latmax];
lonlims = [lonmin lonmax];
figure; worldmap(latlims,lonlims);
set(gcf,'Position',[617, 599, 820, 820])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',18);
land = shaperead('landareas', 'UseGeoCoords', true);
pcolorm(SOCATv2021_grid.latitude(:,:,1),SOCATv2021_grid.longitude(:,:,1),...
    log10(SOCATv2021_grid.p_corr));
geoshow(land, 'FaceColor','w');
c=colorbar; caxis([2.*log10(0.05) 0]); colormap([rgb('darkgreen');rgb('gray')]);
c.Ticks = [2.*log10(0.05) log10(0.05) 0];
c.TickLabels = {'' '0.05' ''};
c.Label.String = 'p-value of correlation (MLD vs. pCO_{2})';
c.Label.FontSize = 18;
exportgraphics(gcf,'/Users/sharp/Desktop/mld_vs_pco2_corr_p.jpg');