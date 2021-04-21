%% Calculate CO2 Flux for individual moorings
% importERA5hourly
% load('/Volumes/2TB Hard Drive/SATELLITE_DATA/ERA5/ERA5_3hourly.mat');

for k = 1:numel(moornames)
    
    % Calculate delta pCO2
    MOORING.(moornames{k}).delta_pCO2 = ...
        MOORING.(moornames{k}).pCO2SW - MOORING.(moornames{k}).pCO2Air;
    % Calculate temporal mean delta pCO2
    MOORING.(moornames{k}).delta_pCO2_PCO2const = ...
        mean(MOORING.(moornames{k}).pCO2SW - MOORING.(moornames{k}).pCO2Air,1,'omitnan');
    MOORING.(moornames{k}).delta_pCO2_PCO2const = repmat(MOORING.(moornames{k}).delta_pCO2_PCO2const,size(MOORING.(moornames{k}).delta_pCO2,1),1);
    MOORING.(moornames{k}).delta_pCO2_PCO2const(isnan(MOORING.(moornames{k}).delta_pCO2)) = NaN;
    
    % Match to wind speed  
    lonidx = ...
        find(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lon,1,'omitnan')) == ...
        min(min(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lon,1,'omitnan')))));
    latidx = ...
        find(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lat,1,'omitnan')) == ...
        min(min(abs(ERA5h.lat-mean(MOORING.(moornames{k}).lat,1,'omitnan')))));
    [~,datidx] = ...
        min(abs(ERA5h.date - datenum(MOORING.(moornames{k}).date)'),[],1);
    MOORING.(moornames{k}).U = squeeze(ERA5h.speed(lonidx,latidx,datidx));
    
    % Calculate gas transfer velocity
    MOORING.(moornames{k}).kw_ERA5 = kgas(MOORING.(moornames{k}).U,660,'W14'); % m/s
    % Calculate gas transfer velocity (U^2 held constant)
    Uconst = sqrt(mean(MOORING.(moornames{k}).U.^2,1,'omitnan'));
    MOORING.(moornames{k}).kw_ERA5_Uconst = kgas(Uconst,660,'W14'); % m/s
    MOORING.(moornames{k}).kw_ERA5_Uconst = repmat(MOORING.(moornames{k}).kw_ERA5_Uconst,size(MOORING.(moornames{k}).kw_ERA5,1),1);

    % Determine K0 (Weiss, R. F., Marine Chemistry 2:203-215, 1974)
    TempK100  = (MOORING.(moornames{k}).sst+273.15)./100;
    lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + MOORING.(moornames{k}).sal .*...
        (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
    MOORING.(moornames{k}).k0_ERA5 = exp(lnK0); % mol/kg-SW/atm

    % Calculate density
    MOORING.(moornames{k}).SA      = gsw_SA_from_SP(MOORING.(moornames{k}).sal,5,MOORING.(moornames{k}).lon,MOORING.(moornames{k}).lat); % Calculate absolute salinity
    MOORING.(moornames{k}).CT      = gsw_CT_from_t(MOORING.(moornames{k}).SA,MOORING.(moornames{k}).sst,5); % Calculate conservative temp
    MOORING.(moornames{k}).DENSITY = gsw_rho(MOORING.(moornames{k}).SA,MOORING.(moornames{k}).CT,5); % Calculate in situ density

    % Calculate flux
    MOORING.(moornames{k}).Fco2_ERA5 = MOORING.(moornames{k}).kw_ERA5 .* MOORING.(moornames{k}).k0_ERA5 .* (MOORING.(moornames{k}).delta_pCO2.*1e-6); % m/s
    MOORING.(moornames{k}).Fco2_ERA5 = MOORING.(moornames{k}).Fco2_ERA5.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
    MOORING.(moornames{k}).Fco2_ERA5 = MOORING.(moornames{k}).Fco2_ERA5.*MOORING.(moornames{k}).DENSITY; % mol/(m^2*yr)
    
    % Calculate flux (U^2 held constant)
    MOORING.(moornames{k}).Fco2_ERA5_Uconst = MOORING.(moornames{k}).kw_ERA5_Uconst .* MOORING.(moornames{k}).k0_ERA5 .* (MOORING.(moornames{k}).delta_pCO2.*1e-6); % m/s
    MOORING.(moornames{k}).Fco2_ERA5_Uconst = MOORING.(moornames{k}).Fco2_ERA5_Uconst.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
    MOORING.(moornames{k}).Fco2_ERA5_Uconst = MOORING.(moornames{k}).Fco2_ERA5_Uconst.*MOORING.(moornames{k}).DENSITY; % mol/(m^2*yr)

    % Calculate flux (delta pCO2 held constant)
    MOORING.(moornames{k}).Fco2_ERA5_PCO2const = MOORING.(moornames{k}).kw_ERA5 .* MOORING.(moornames{k}).k0_ERA5 .* (MOORING.(moornames{k}).delta_pCO2_PCO2const.*1e-6); % m/s
    MOORING.(moornames{k}).Fco2_ERA5_PCO2const = MOORING.(moornames{k}).Fco2_ERA5_PCO2const.*(60.*60.*24.*365.25); % (m*mol)/(yr*kg)
    MOORING.(moornames{k}).Fco2_ERA5_PCO2const = MOORING.(moornames{k}).Fco2_ERA5_PCO2const.*MOORING.(moornames{k}).DENSITY; % mol/(m^2*yr)

    % Determine monthly means
    monthly_dates = datenum([[repelem([1998:2019]',12,1);2020] [repmat([1:12]',22,1);1] ones(22*12+1,1)]);
    [count,bins,num] = histcounts(datenum(MOORING.(moornames{k}).date),monthly_dates);
    MOORING.(moornames{k}).Fco2_ERA5_mon = accumarray(num,MOORING.(moornames{k}).Fco2_ERA5,[264 1],@mean,NaN);
    MOORING.(moornames{k}).Fco2_ERA5_Uconst_mon = accumarray(num,MOORING.(moornames{k}).Fco2_ERA5_Uconst,[264 1],@mean,NaN);
    MOORING.(moornames{k}).Fco2_ERA5_PCO2const_mon = accumarray(num,MOORING.(moornames{k}).Fco2_ERA5_PCO2const,[264 1],@mean,NaN);
    
    % Determine RMSEs for each component
    MOORING.(moornames{k}).RMSE_ERA5_Uconst = ...
        sqrt(mean((MOORING.(moornames{k}).Fco2_ERA5_Uconst - MOORING.(moornames{k}).Fco2_ERA5).^2,1,'omitnan'));
    MOORING.(moornames{k}).RMSE_ERA5_PCO2const = ...
        sqrt(mean((MOORING.(moornames{k}).Fco2_ERA5_PCO2const - MOORING.(moornames{k}).Fco2_ERA5).^2,1,'omitnan'));

    % Determine RMSE ratio
    MOORING.(moornames{k}).RMSE_ratio_ERA5 = MOORING.(moornames{k}).RMSE_ERA5_Uconst./MOORING.(moornames{k}).RMSE_ERA5_PCO2const;

    % Plot monthly means
    figure; hold on;
    set(gcf,'Position',[100, 400, 1680, 420])
    set(gca,'fontsize',16)
    title((moornames{k}),'fontsize',18);
    p1=plot(monthly_dates(1:end-1),MOORING.(moornames{k}).Fco2_ERA5_mon,'k','linewidth',3);
    p2=plot(monthly_dates(1:end-1),MOORING.(moornames{k}).Fco2_ERA5_Uconst_mon,'b','linewidth',3);
    p3=plot(monthly_dates(1:end-1),MOORING.(moornames{k}).Fco2_ERA5_PCO2const_mon,'r','linewidth',3);
    text(min(monthly_dates(1:end-1)),max(MOORING.(moornames{k}).Fco2_ERA5_mon),strcat('|U|^{2} effect / \DeltapCO_{2} effect = ',num2str(MOORING.(moornames{k}).RMSE_ratio_ERA5)),'fontsize',16)
    xlabel('Year');
    ylabel('CO2 Flux');
    datetick('x','yyyy');
    legend([p1 p2 p3],{'CO2 flux'...
                   'U^{2} held constant'...
                   '\DeltapCO2 held constant'},'location','southeast');
    exportgraphics(gcf,strcat('/Users/sharp/Desktop/Monthly_Wind_Components_',moornames{k},'.jpg'));

    
    % Plot individual data points
    figure; hold on;
    set(gcf,'Position',[100, 400, 1680, 420])
    set(gca,'fontsize',16)
    title((moornames{k}),'fontsize',18);
    p1=plot(datenum(MOORING.(moornames{k}).date),MOORING.(moornames{k}).Fco2_ERA5,'k.','linewidth',3);
    p2=plot(datenum(MOORING.(moornames{k}).date),MOORING.(moornames{k}).Fco2_ERA5_Uconst,'b.','linewidth',3);
    p3=plot(datenum(MOORING.(moornames{k}).date),MOORING.(moornames{k}).Fco2_ERA5_PCO2const,'r.','linewidth',3);
    text(min(datenum(MOORING.(moornames{k}).date)),max(MOORING.(moornames{k}).Fco2_ERA5),strcat('|U|^{2} effect / \DeltapCO_{2} effect = ',num2str(MOORING.(moornames{k}).RMSE_ratio_ERA5)),'fontsize',16)
    datetick('x','yyyy');
    legend([p1 p2 p3],{'CO2 flux'...
                   'U^{2} held constant'...
                   '\DeltapCO2 held constant'},'location','southeast');
    exportgraphics(gcf,strcat('/Users/sharp/Desktop/Wind_Components_',moornames{k},'.jpg'));

end