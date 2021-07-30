%% Calculate CO2 Flux for individual moorings
% importERA5hourly
% load('/Volumes/2TB Hard Drive/SATELLITE_DATA/ERA5/ERA5_3hourly.mat');

moornames  = {'WA' 'CCE1' 'CCE2' 'NH10' 'LaPush'};
moornames2 = {'Cape Elizabeth' 'CCE1' 'CCE2' 'NH10' 'Chá bâ'};

meanflux = nan(length(moornames),4);
cumflux = nan(length(moornames),4);

for k = 3 % CCE2 only 1:numel(moornames)
    
    yr = 2015;
    
    % Index to desired timespan
    dateidx = ...
        MOORING.(moornames{k}).date(:,1) == yr;

    % extract datenum from mooring observations
    [h,m,s] = hms(MOORING.(moornames{k}).time(dateidx));
    date_h = datenum([MOORING.(moornames{k}).date(dateidx,1:3),h,m,s]);   
    
    % Calculate delta pCO2
    MOORING.(moornames{k}).delta_pCO2 = ...
        MOORING.(moornames{k}).pCO2SW(dateidx) - MOORING.(moornames{k}).pCO2Air(dateidx); % uatm

    % randomly select one delta pco2 value per month
    MOORING.(moornames{k}).delta_pCO2_rnd = nan(12,1);
    rng(2); %sets random observations that are selected
    for n = 1:12
        pco2_subset = MOORING.(moornames{k}).delta_pCO2(MOORING.(moornames{k}).date(dateidx,2)==n);
        choice = NaN;
        while isnan(choice)
            idx = randi(length(pco2_subset));
            choice = pco2_subset(idx);
        end
        MOORING.(moornames{k}).delta_pCO2_rnd(n) = choice;
    end
    
    % Match to wind speed  
    lonidx = ...
        find(abs(ERA5h.lon-mean(MOORING.(moornames{k}).lon(dateidx),1,'omitnan')) == ...
        min(min(abs(ERA5h.lon-mean(MOORING.(moornames{k}).lon(dateidx),1,'omitnan')))));
    latidx = ...
        find(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lat(dateidx),1,'omitnan')) == ...
        min(min(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lat(dateidx),1,'omitnan')))));
    [~,datidx] = ...
        min(abs(ERA5h.date - datenum(MOORING.(moornames{k}).date(dateidx))'),[],1);
    MOORING.(moornames{k}).U = squeeze(ERA5h.speed(lonidx,latidx,datidx));

    % Calculate gas transfer velocity for mooring with kgas function
    MOORING.(moornames{k}).sch = CO2flux_Schmidt_W14(MOORING.(moornames{k}).sst(dateidx),'CO2');
    MOORING.(moornames{k}).kw_ERA5 = kgas(MOORING.(moornames{k}).U,MOORING.(moornames{k}).sch,'Ho06'); % m/s
    MOORING.(moornames{k}).kw_ERA5 = MOORING.(moornames{k}).kw_ERA5.*(60.*60.*24.*365.25./12); % m/h
    MOORING.(moornames{k}).kw_ERA5_hourly = MOORING.(moornames{k}).kw_ERA5./(24.*365.25./12); % m/h

    % Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
    TempK100  = (MOORING.(moornames{k}).sst(dateidx)+273.15)./100;
    lnK0 = -58.0931 + 90.5069 ./ TempK100 + 22.2940 .* log(TempK100) + MOORING.(moornames{k}).sal(dateidx) .*...
        (0.027766 - 0.025888 .* TempK100 + 0.0050578 .* TempK100 .^2);
    MOORING.(moornames{k}).k0_ERA5 = exp(lnK0).*1e3; % mmol/(l*atm)
    MOORING.(moornames{k}).k0_ERA5 = MOORING.(moornames{k}).k0_ERA5.*1000; % mmol/(m^3 * atm)
    
    % Match to wind speed  
    lonidx = ...
        find(abs(SOCATv2021_grid.lon(:,1,1)-mean(MOORING.(moornames{k}).lon(dateidx),1,'omitnan')) == ...
        min(min(abs(SOCATv2021_grid.lon(:,1,1)-mean(MOORING.(moornames{k}).lon(dateidx),1,'omitnan')))));
    latidx = ...
        find(abs(SOCATv2021_grid.lat(1,:,1)'-mean(MOORING.(moornames{k}).lat(dateidx),1,'omitnan')) == ...
        min(min(abs(SOCATv2021_grid.lat(1,:,1)'-mean(MOORING.(moornames{k}).lat(dateidx),1,'omitnan'))))); 
    
    % Calculate flux
    MOORING.(moornames{k}).Fco2_ERA5_rnd = ...
        squeeze(SOCATv2021_grid.kw_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))) .* ...
        squeeze(SOCATv2021_grid.k0_ERA5(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))) .* ...
        (MOORING.(moornames{k}).delta_pCO2_rnd.*1e-6);  % mmol C/(m^2 * h)
    MOORING.(moornames{k}).Fco2_ERA5_hourly = ...
        MOORING.(moornames{k}).kw_ERA5_hourly .* ...
        MOORING.(moornames{k}).k0_ERA5 .* ...
        (MOORING.(moornames{k}).delta_pCO2.*1e-6);  % mmol C/(m^2 * h)

    % Determine monthly means
    monthly_dates = datenum([[repelem([1998:2019]',12,1);2020] [repmat([1:12]',22,1);1] ones(22*12+1,1)]);
    [count,bins,num] = histcounts(datenum(MOORING.(moornames{k}).date(dateidx,:)),monthly_dates);
    MOORING.(moornames{k}).Fco2_ERA5_mon = accumarray(num,MOORING.(moornames{k}).Fco2_ERA5,[276 1],@nanmean,NaN); % mmol C / m^2 hr
    
    % calculate dates
    date = datenum([repmat(1998,size(SOCATv2021_grid.month_since_1998))...
                    SOCATv2021_grid.month_since_1998...
                    repmat(15,size(SOCATv2021_grid.month_since_1998))]);

    % initialize figure
    pos = [0 0 1 0.5];
    titlesz = 28;
    fontsz = 22;
    figure; box on; hold on;
    set(gcf,'units','normalized','outerposition',pos);
    set(gca,'fontsize',fontsz)
    title((moornames2{k}),'fontsize',titlesz);
    
    % Plot monthly means
    s1=scatter(date_h,MOORING.(moornames{k}).Fco2_ERA5_hourly,...
        'MarkerFaceColor',rgb('grey'),'MarkerEdgeColor',rgb('grey'));
    p2=plot(date((yr-1998)*12+1:(yr-1998)*12+12),MOORING.(moornames{k}).Fco2_ERA5_rnd,...
        'k','linewidth',4);
    s2=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),MOORING.(moornames{k}).Fco2_ERA5_rnd,100,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');  % mmol C / m^2 hr
    p3=plot(date((yr-1998)*12+1:(yr-1998)*12+12),...
        squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),...
        '-','linewidth',4,'Color',[0 0 0.9]); % mmol C / m^2 hr
    s3=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),...
        squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))),100,...
        'MarkerFaceColor',[0 0 0.9],'MarkerEdgeColor','k'); % mmol C / m^2 hr
%     p4=plot(date((yr-1998)*12+1:(yr-1998)*12+12),squeeze(LAR.Fco2_ERA5_hourly(lonidx,latidx,((yr-1998)*12+1:(yr-1998)*12+12))),...
%         '-','Color',[0 0.7 0],'linewidth',4); % mmol C / m^2 hr
%     s4=scatter(date((yr-1998)*12+1:(yr-1998)*12+12),squeeze(LAR.Fco2_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),100,...
%         'MarkerFaceColor',[0 0.7 0],'MarkerEdgeColor','k'); % mmol C / m^2 hr
    plot([735960 736330],[0 0],'k--','linewidth',1)
    ylabel('{\itF}_{CO2} (mmol C m^{-2} hr^{-1})','fontsize',fontsz);
    ylim([-0.4 0.4]);
    %datetick('x','m','keeplimits');
    xticks(date((yr-1998)*12+1:(yr-1998)*12+12));
    xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
    xlim([min(datenum(MOORING.(moornames{k}).date(dateidx,:))) ...
        max(datenum(MOORING.(moornames{k}).date(dateidx,:)))]);
    legend([s1 p2 p3],{'3-hourly flux (Mooring)'...
                   'random monthly flux (Mooring)'...
                   'monthly flux (RFR-CCS)'},'location','northeast');
    exportgraphics(gcf,strcat('/Users/sharp/Desktop/Flux_',moornames{k},'_2.jpg'));
    
%     meanflux(k,1) = ...
%         mean(MOORING.(moornames{k}).Fco2_ERA5,'omitnan').*24.*365.25/12./1e3 % mol C / m^2 yr
% %     index = date > min(datenum(MOORING.(moornames{k}).date)) & ...
% %             date < max(datenum(MOORING.(moornames{k}).date));
%     meanflux(k,2) = ...
%         mean(MOORING.(moornames{k}).Fco2_ERA5_mon((yr-1998)*12+1:(yr-1998)*12+12),'omitnan').*24.*365.25/12./1e3 % mol C / m^2 yr
%     meanflux(k,3) = ...
%         mean(squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),'omitnan').*24.*365.25/12./1e3 % mol C / m^2 yr
%     meanflux(k,4) = ...
%         mean(squeeze(LAR.Fco2_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12)),'omitnan').*24.*365.25/12./1e3 % mol C / m^2 yr
%     
    % interpolate over missing 3-hourly observations
    MOORING.(moornames{k}).Fco2_ERA5_interp = ...
        interp1(date_h(~isnan(MOORING.(moornames{k}).Fco2_ERA5)),...
        MOORING.(moornames{k}).Fco2_ERA5(~isnan(MOORING.(moornames{k}).Fco2_ERA5)),date_h);
    
        sum(MOORING.(moornames{k}).Fco2_ERA5_interp,'omitnan').*3 % mol C / m^2 yr
        sum(MOORING.(moornames{k}).Fco2_ERA5_rnd).*(3.*2918/12) % mol C / m^2 yr
        sum(squeeze(SOCATv2021_grid.Fco2_RF_ERA5_hourly(lonidx,latidx,(yr-1998)*12+1:(yr-1998)*12+12))).*(3.*2918/12) % mol C / m^2 yr


end

%% Plot bar graph of amplitudes
% figure; hold on; box on;
% set(gca,'fontsize',18);
% set(gcf,'Position',[617, 599, 800, 400])
% x_ax = categorical(moornames2);
% bargra = bar(x_ax,meanflux,'FaceColor','flat');
% bargra(1).CData = rgb('grey');
% bargra(2).CData = rgb('black');
% bargra(3).CData = [0 0 0.9];
% bargra(4).CData = [0 0.7 0];
% ylabel('Mean flux (\mumol C m^{-2} hr^{-1})','fontsize',18);
% ylim([-8 2]);
% legend({'mooring data (3-hourly)' 'mooring data (monthly)' 'RFR-CCS (monthly)' 'L17 (monthly)'},...
%     'location','southwest','fontsize',14)
% exportgraphics(gcf,strcat('/Users/sharp/Desktop/Mooring_Flux_Bar.jpg'));

