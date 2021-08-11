%% Calculate CO2 Flux for the CCE2 mooring
load('Data/CCE2_winds_hourly.mat');

% Define year of interest
yr = 2015;

% Index to desired timespan
dateidx = MOORING.CCE2.date(:,1) == yr;

% extract datenum from mooring observations
[h,m,s] = hms(MOORING.CCE2.time(dateidx));
date_h = datenum([MOORING.CCE2.date(dateidx,1:3),h,m,s]);   

% Calculate delta pCO2
MOORING.CCE2.delta_pCO2 = ...
    MOORING.CCE2.pCO2SW(dateidx) - MOORING.CCE2.pCO2Air(dateidx); % uatm

% Match RFR-CCS lat/lon to mooring lat/lon
lonidx = ...
    find(abs(SOCATv2021_grid.lon(:,1,1)-mean(MOORING.CCE2.lon(dateidx),1,'omitnan')) == ...
    min(min(abs(SOCATv2021_grid.lon(:,1,1)-mean(MOORING.CCE2.lon(dateidx),1,'omitnan')))));
latidx = ...
    find(abs(SOCATv2021_grid.lat(1,:,1)'-mean(MOORING.CCE2.lat(dateidx),1,'omitnan')) == ...
    min(min(abs(SOCATv2021_grid.lat(1,:,1)'-mean(MOORING.CCE2.lat(dateidx),1,'omitnan'))))); 

MOORING.CCE2.delta_pCO2_rnd = nan(1000,12);
MOORING.CCE2.Fco2_ERA5_rnd = nan(1000,12);
for t=1:100000
% randomly select one delta pco2 value per month
%rng(2); %sets random observations that are selected
for n = 1:12
    pco2_subset = MOORING.CCE2.delta_pCO2(MOORING.CCE2.date(dateidx,2)==n);
    choice = NaN;
    while isnan(choice)
        idx = randi(length(pco2_subset));
        choice = pco2_subset(idx);
    end
    MOORING.CCE2.delta_pCO2_rnd(t,n) = choice;
end
MOORING.CCE2.Fco2_ERA5_rnd(t,:) = ...
    squeeze(SOCATv2021_grid.kw_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12)))' .* ...
    squeeze(SOCATv2021_grid.k0_ERA5(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12)))' .* ...
    (MOORING.CCE2.delta_pCO2_rnd(t,:).*1e-6);  % mmol C/(m^2 * h)
end
MOORING.CCE2.Fco2_ERA5_rnd_mean = mean(MOORING.CCE2.Fco2_ERA5_rnd,1);
MOORING.CCE2.Fco2_ERA5_rnd_std = std(MOORING.CCE2.Fco2_ERA5_rnd,1);

% Match to time of wind speed
datidx = nan(length(date_h),1);
for n = 1:length(date_h)
datidx(n) = find(abs(CCE2_winds_hourly_date - date_h(n)) == ...
    min(abs(CCE2_winds_hourly_date - date_h(n)),[],1));
end

CCE2_winds_hourly = CCE2_winds_hourly(datidx);

% Calculate gas transfer velocity for mooring with kgas function
MOORING.CCE2.sch = CO2flux_Schmidt_W14(MOORING.CCE2.sst(dateidx),'CO2');
MOORING.CCE2.kw_ERA5 = kgas(CCE2_winds_hourly,MOORING.CCE2.sch,'Ho06'); % m/s
MOORING.CCE2.kw_ERA5 = MOORING.CCE2.kw_ERA5.*(60.*60.*24.*365.25./12); % m/h
MOORING.CCE2.kw_ERA5_hourly = MOORING.CCE2.kw_ERA5./(24.*365.25./12); % m/h

% Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
TempK100  = (MOORING.CCE2.sst(dateidx)+273.15)./100;
lnK0 = -58.0931 + 90.5069 ./ TempK100 + 22.2940 .* log(TempK100) + MOORING.CCE2.sal(dateidx) .*...
    (0.027766 - 0.025888 .* TempK100 + 0.0050578 .* TempK100 .^2);
MOORING.CCE2.k0_ERA5 = exp(lnK0).*1e3; % mmol/(l*atm)
MOORING.CCE2.k0_ERA5 = MOORING.CCE2.k0_ERA5.*1000; % mmol/(m^3 * atm)

% Calculate flux
% MOORING.CCE2.Fco2_ERA5_rnd = ...
%     squeeze(SOCATv2021_grid.kw_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))) .* ...
%     squeeze(SOCATv2021_grid.k0_ERA5(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))) .* ...
%     (MOORING.CCE2.delta_pCO2_rnd.*1e-6);  % mmol C/(m^2 * h)
MOORING.CCE2.Fco2_ERA5_hourly = ...
    MOORING.CCE2.kw_ERA5_hourly .* ...
    MOORING.CCE2.k0_ERA5 .* ...
    (MOORING.CCE2.delta_pCO2.*1e-6);  % mmol C/(m^2 * h)
% interpolate over missing 3-hourly observations
MOORING.CCE2.Fco2_ERA5_interp = ...
    interp1(date_h(~isnan(MOORING.CCE2.Fco2_ERA5_hourly)),...
    MOORING.CCE2.Fco2_ERA5_hourly(~isnan(MOORING.CCE2.Fco2_ERA5_hourly)),date_h);

% Determine monthly means
monthly_dates = datenum([[repelem([1998:2019]',12,1);2020] [repmat([1:12]',22,1);1] ones(22*12+1,1)]);
[count,bins,num] = histcounts(datenum(MOORING.CCE2.date(dateidx,:)),monthly_dates);
MOORING.CCE2.Fco2_ERA5_mon = accumarray(num,MOORING.CCE2.Fco2_ERA5_interp,[276 1],@nanmean,NaN); % mmol C / m^2 hr

% Determine cumulative annual flux
MOORING.CCE2.Fco2_ERA5_hourly_cum = ...
    sum(MOORING.CCE2.Fco2_ERA5_interp,'omitnan').*3; % mmol C / m^2 yr
SOCATv2021_grid.Fco2_RF_ERA5_hourly_cum = ...
    sum(squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,...
    (yr-1998)*12+1:(yr-1998)*12+12))).*(3.*2918/12); % mmol C / m^2 yr
MOORING.CCE2.Fco2_ERA5_rnd_mean_cum = ...
    sum(MOORING.CCE2.Fco2_ERA5_rnd_mean).*(3.*2918/12); % mmol C / m^2 yr
MOORING.CCE2.Fco2_ERA5_rnd_std_cum = ...
    sum(MOORING.CCE2.Fco2_ERA5_rnd_std).*(3.*2918/12); % mmol C / m^2 yr

% calculate dates
date = datenum([repmat(1998,size(SOCATv2021_grid.month_since_1998))...
                SOCATv2021_grid.month_since_1998...
                repmat(15,size(SOCATv2021_grid.month_since_1998))]);

%% create plot
% initialize figure
pos = [0 0 1 0.5];
titlesz = 28;
fontsz = 22;
figure;
set(gcf,'units','normalized','outerposition',pos);
ax1=axes('position',[0.05 0.1 0.73 0.8],'units','normalized','box','on');
set(ax1,'fontsize',fontsz)
title('CCE2','fontsize',titlesz);
hold on;

% Plot monthly means
s1=scatter(date_h,MOORING.CCE2.Fco2_ERA5_hourly,...
    'MarkerFaceColor',rgb('grey'),'MarkerEdgeColor',rgb('grey'));
p2=fill([date((yr-1998)*12+1:(yr-1998)*12+12);flipud(date((yr-1998)*12+1:(yr-1998)*12+12))],...
    [MOORING.CCE2.Fco2_ERA5_rnd_mean+MOORING.CCE2.Fco2_ERA5_rnd_std ...
    fliplr(MOORING.CCE2.Fco2_ERA5_rnd_mean-MOORING.CCE2.Fco2_ERA5_rnd_std)]',...
    [0 0.7 0],'FaceAlpha',0.5,'LineStyle','none');
% p2=plot(date((yr-1998)*12+1:(yr-1998)*12+12),MOORING.CCE2.Fco2_ERA5_rnd_mean+MOORING.CCE2.Fco2_ERA5_rnd_std,...
%     'k','linewidth',4);
% p22=plot(date((yr-1998)*12+1:(yr-1998)*12+12),MOORING.CCE2.Fco2_ERA5_rnd_mean-MOORING.CCE2.Fco2_ERA5_rnd_std,...
%     'k','linewidth',4);
% s2=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),MOORING.CCE2.Fco2_ERA5_rnd_mean,100,...
%     'MarkerFaceColor','k','MarkerEdgeColor','k');  % mmol C / m^2 hr
p3=plot(date((yr-1998)*12+1:(yr-1998)*12+12),...
    squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),...
    '-','linewidth',6,'Color',[0 0 0.9]); % mmol C / m^2 hr
s3=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),...
    squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))),200,...
    'MarkerFaceColor',[0 0 0.9],'MarkerEdgeColor','k'); % mmol C / m^2 hr
%     p4=plot(date((yr-1998)*12+1:(yr-1998)*12+12),squeeze(LAR.Fco2_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))),...
%         '-','Color',[0 0.7 0],'linewidth',4); % mmol C / m^2 hr
%     s4=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),squeeze(LAR.Fco2_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),100,...
%         'MarkerFaceColor',[0 0.7 0],'MarkerEdgeColor','k'); % mmol C / m^2 hr
plot([735960 736330],[0 0],'k--','linewidth',1);
ylabel('{\itF}_{CO2} (mmol C m^{-2} hr^{-1})','fontsize',fontsz);
%ylim([-0.4 0.4]);
%datetick('x','m','keeplimits');
xticks(date((yr-1998)*12+1:(yr-1998)*12+12));
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
xlim([min(datenum(MOORING.CCE2.date(dateidx,:))) ...
    max(datenum(MOORING.CCE2.date(dateidx,:)))]);
legend([s1 p2 p3],{'3-hourly flux (Mooring)'...
               'random monthly flux (Mooring)'...
               'monthly flux (RFR-CCS)'},'location','northeast');

% Plot cumulative annual flux
ax2 = axes('position',[0.84 0.15 0.15 0.8],'units','normalized','box','on');
x_ax = categorical({'3-hourly' 'random' 'RFR-CCS'});
[h,e] = barwitherr([NaN NaN MOORING.CCE2.Fco2_ERA5_rnd_std_cum],...
    [MOORING.CCE2.Fco2_ERA5_hourly_cum ...
    SOCATv2021_grid.Fco2_RF_ERA5_hourly_cum ...
    MOORING.CCE2.Fco2_ERA5_rnd_mean_cum]);
% bargra = bar(x_ax,[MOORING.CCE2.Fco2_ERA5_hourly_cum ...
%     MOORING.CCE2.Fco2_ERA5_rnd_mean_cum ...
%     SOCATv2021_grid.Fco2_RF_ERA5_hourly_cum],'FaceColor','flat');
h.CData = [rgb('grey');[0 0 0.9];[0 0 0.9]]
h(2).CData = rgb('blue')
set(h,'FaceColor',[rgb('grey');[0 0 0.9];[0 0 0.9]]);
set(h(2),'FaceColor',[0 0 0.9]);
set(h(3),'FaceColor',[0 0.7 0]);
set(ax2,'fontsize',fontsz);
ylabel('Cumulative annual flux (mmol C m^{-2} yr^{-1})','fontsize',fontsz);
ylim([-800 0]);

exportgraphics(gcf,strcat('/Users/sharp/Desktop/Flux_CCE2.jpg'));

