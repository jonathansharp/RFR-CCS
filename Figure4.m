%% This script produces Figure 4 from Sharp et al. (in prep)

% Import SOCATv4 gridded files
importSOCATv4gridded
importSOCATv4gridded_coastal

SOCATv4_gridded.distance_from_shore = ...
    dist2coast(SOCATv4_gridded.lat,SOCATv4_gridded.lon);
SOCATv4_gridded.distance_from_shore = ...
    repmat(SOCATv4_gridded.distance_from_shore,1,1,...
    size(SOCATv4_gridded.month_since_1998,1));

%% Initialize figure

labelsz = 24;
fontsz = 20;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.75 0.75])
tl = axes('Position',[0.06 0.59 0.37 0.38],'Box','on');
tr = axes('Position',[0.53 0.59 0.37 0.38],'Box','on');
bl = axes('Position',[0.06 0.09 0.37 0.38],'Box','on');
br = axes('Position',[0.53 0.09 0.37 0.38],'Box','on');

%% Interpolate values to open-ocean grid cells
for m = 1:12
    % RFR climatology
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,SOCATv2021_grid.pco2_RF_clim_2015(:,:,m));
    SOCATv2021_grid.pco2_RF_clim_open_2015(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    % Landschutzer product
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,LAND.pCO2(:,:,m));
    LAND.pCO2_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
    
end
for m = 1:264
    % RFR monthly product
    interp = griddedInterpolant(SOCATv2021_grid.lon,SOCATv2021_grid.lat,SOCATv2021_grid.pco2_RF(:,:,m));
    SOCATv2021_grid.pco2_RF_open(:,:,m) = ...
        interp(360+SOCATv4_gridded.lon,SOCATv4_gridded.lat);
end

%% Bin values for 2d histograms
% RFR values vs SOCAT grid cells (climatological)
[counts_1,log_dist_1,x_bin_centers_1,y_bin_centers_1,Pcorr_1,total_1] = ...
    bin_values(SOCATv4_gridded.pCO2_mon_mean_clim,...
    SOCATv4_gridded_coastal.pCO2_mon_mean_clim,...
    SOCATv2021_grid.pco2_RF_clim_open_2015,...
    SOCATv2021_grid.pco2_RF_clim_2015,...
    SOCATv4_gridded.distance_from_shore,...
    SOCATv2021_grid.distance_from_shore(:,:,1:12),...
    nan,1);
% Landschutzer values vs SOCAT grid cells (climatological)
[counts_2,log_dist_2,x_bin_centers_2,y_bin_centers_2,Pcorr_2,total_2] = ...
    bin_values(SOCATv4_gridded.pCO2_mon_mean_clim,...
    SOCATv4_gridded_coastal.pCO2_mon_mean_clim,...
    LAND.pCO2_open,...
    LAND.pCO2,...
    SOCATv4_gridded.distance_from_shore,...
    SOCATv2021_grid.distance_from_shore(:,:,1:12),...
    nan,2);
% RFR values vs SOCAT grid cells (monthly coastal)
[counts_3,log_dist_3,x_bin_centers_3,y_bin_centers_3,Pcorr_3,total_3] = ...
    bin_values([],...
    SOCATv4_gridded_coastal.pCO2_mon_mean,...
    [],...
    SOCATv2021_grid.pco2_RF(:,:,1:216),...
    [],...
    SOCATv2021_grid.distance_from_shore(:,:,1:216),...
    LAR.pCO2,3);
% Laruelle values vs SOCAT grid cells (monthly coastal)
[counts_4,log_dist_4,x_bin_centers_4,y_bin_centers_4,Pcorr_4,total_4] = ...
    bin_values([],...
    SOCATv4_gridded_coastal.pCO2_mon_mean,...
    [],...
    LAR.pCO2,...
    [],...
    SOCATv2021_grid.distance_from_shore(:,:,1:216),...
    nan,4);

%% Plot RFR values vs SOCAT grid cells (climatological)
axes(tl); hold on;
h=image('XData',x_bin_centers_1,'YData',y_bin_centers_1,'CData',log_dist_1');
h.CDataMapping = 'scaled';
log_counts_1 = log10(counts_1');
log_counts_1(log_counts_1==-Inf) = 0;
set(h,'AlphaData',log_counts_1);
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([150 650]);
ylabel('RFR-CCS-clim {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([150 600]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
% xticks([200 400 600 800]);
% yticks([200 400 600 800]);
text(200,550,strcat('R^{2} =',{' '},...
    num2str(round(Pcorr_1,2))),'fontsize',labelsz);
text(210,500,strcat('{\itn} =',{' '},...
    num2str(total_1)),'fontsize',labelsz-6);
total_1

%% Plot Landschutzer values vs SOCAT grid cells
axes(tr); hold on;
h=image('XData',x_bin_centers_2,'YData',y_bin_centers_2,'CData',log_dist_2);
h.CDataMapping = 'scaled';
log_counts_2 = log10(counts_2');
log_counts_2(log_counts_2==-Inf) = 0;
set(h,'AlphaData',log_counts_2);
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly mean {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([150 650]);
ylabel('L20 {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([150 600]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
% xticks([200 400 600 800]);
% yticks([200 400 600 800]);
text(200,550,strcat('R^{2} =',{' '},...
    num2str(round(Pcorr_2,2))),'fontsize',labelsz);
text(210,500,strcat('{\itn} =',{' '},...
    num2str(total_2)),'fontsize',labelsz-6);
total_2

%% Plot RFR values vs SOCAT grid cells (monthly coastal)
axes(bl); hold on;
h=image('XData',x_bin_centers_3,'YData',y_bin_centers_3,'CData',log_dist_3);
h.CDataMapping = 'scaled';
log_counts_3 = log10(counts_3');
log_counts_3(log_counts_3==-Inf) = 0;
set(h,'AlphaData',log_counts_3);
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([150 650]);
ylabel('RFR-CCS {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([150 600]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
% xticks([200 400 600 800]);
% yticks([200 400 600 800]);
text(200,550,strcat('R^{2} =',{' '},...
    num2str(round(Pcorr_3,2))),'fontsize',labelsz);
text(210,500,strcat('{\itn} =',{' '},...
    num2str(total_3)),'fontsize',labelsz-6);
total_3

%% Plot Laruelle values vs SOCAT grid cells
axes(br); hold on;
h=image('XData',x_bin_centers_4,'YData',y_bin_centers_4,'CData',log_dist_4);
h.CDataMapping = 'scaled';
log_counts_4 = log10(counts_4');
log_counts_4(log_counts_4==-Inf) = 0;
set(h,'AlphaData',log_counts_4);
set(gca,'fontsize',fontsz);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('SOCATv4 gridded monthly {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
xlim([150 650]);
ylabel('L17 {\itp}CO_{2(sw)} (\muatm)','fontsize',labelsz);
ylim([150 600]);
box('on');
colormap(cmocean('haline'));
caxis([1 3]);
% xticks([200 400 600 800]);
% yticks([200 400 600 800]);
text(200,550,strcat('R^{2} =',{' '},...
    num2str(round(Pcorr_4,2))),'fontsize',labelsz);
text(210,500,strcat('{\itn} =',{' '},...
    num2str(total_4)),'fontsize',labelsz-6);
total_4

%% Colorbar
ax=axes('Position',[0.06 0.05 0.9 0.9],'Box','off');
set(gca,'fontsize',fontsz);
c=colorbar;
caxis([1 3]);
colormap(cmocean('haline'));
c.Label.String = 'log_{10}(km from shore)';
c.Label.FontSize = labelsz;
ax.Visible = 'off';

annotation('textbox',[0.4, 0.95, 0, 0],'string','a','fontsize',24)
annotation('textbox',[0.86, 0.95, 0, 0],'string','b','fontsize',24)
annotation('textbox',[0.4, 0.45, 0, 0],'string','c','fontsize',24)
annotation('textbox',[0.86, 0.45, 0, 0],'string','d','fontsize',24)

exportgraphics(gcf,'Figures/Figure4.png');

%% embedded function
function [counts,log_dist,x_bin_centers,y_bin_centers,Pcorr,total] = ...
    bin_values(x_open,x_coastal,y_open,y_coastal,dist_open,dist_coastal,extra,num)

bins = -2.5:5:1202.5;
IndexOpen = ~isnan(x_open) & ~isnan(y_open);
IndexCoast = ~isnan(x_coastal) & ~isnan(y_coastal);
if num == 3
IndexCoast = ~isnan(extra) & ~isnan(x_coastal) & ~isnan(y_coastal);
end
[counts,bin_centers] = ...
    hist3([[x_open(IndexOpen);x_coastal(IndexCoast)],...
    [y_open(IndexOpen);y_coastal(IndexCoast)]],...
    'Edges',{bins bins});
[~,~,Xnum] = ...
    histcounts([x_open(IndexOpen);x_coastal(IndexCoast)],bins);
[~,~,Ynum] = ...
    histcounts([y_open(IndexOpen);y_coastal(IndexCoast)],bins);
log_dist = accumarray([Xnum Ynum],...
    [log10(dist_open(IndexOpen));log10(dist_coastal(IndexCoast))],...
    [length(bins) length(bins)],@nanmean);
x_bin_centers = bin_centers{1};
y_bin_centers = bin_centers{2};
Pcorr = corr([x_open(IndexOpen);x_coastal(IndexCoast)],...
    [y_open(IndexOpen);y_coastal(IndexCoast)]);
total = sum(~isnan([x_open(IndexOpen);x_coastal(IndexCoast)]));

end