% Mooring omission analysis (Fig 3)

MOORING.SEAK = MOORING.Southeast;
MOORING.CB06 = MOORING.CB_06;

omitmoors2 = {'WA' 'CCE1' 'CCE2' 'NH10' 'LaPush'};
omitmoor_names = {'Cape Elizabeth' 'CCE1' 'CCE2' 'NH10' 'Châ bá'};

Mean_table = nan(max(size(omitmoors2)),4);
Amplitude_table = nan(max(size(omitmoors2)),4);
R2_table = nan(max(size(omitmoors2)),3);
Bias_table = nan(max(size(omitmoors2)),3);
RMSE_table = nan(max(size(omitmoors2)),3);

%% For each omitted mooring
for k=1:numel(omitmoors2)
    
    % Match latitude in climatology
    latidx = ...
        find(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(omitmoors2{k}).lat)) == ...
        min(min(abs(SOCATv2021_grid.latitude(1,:,1)-nanmean(MOORING.(omitmoors2{k}).lat)))));
    lonidx = ...
        find(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(omitmoors2{k}).lon)) == ...
        min(min(abs(SOCATv2021_grid.longitude(:,1,1)-nanmean(MOORING.(omitmoors2{k}).lon)))));
    
    % Determine mooring month since 1998
    MOORING.(omitmoors2{k}).month_since_1998 = ...
        (MOORING.(omitmoors2{k}).date(:,1)-1998).*12 + MOORING.(omitmoors2{k}).date(:,2);
    
    % Determine means of each month since 1998
    MOORING.(omitmoors2{k}).pCO2SW_mon = nan(max(SOCATv2021_grid.month_since_1998),1);
    for m=1:max(SOCATv2021_grid.month_since_1998)
        MOORING.(omitmoors2{k}).month_since_1998_mon(m) = m;
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(m) = mean(MOORING.(omitmoors2{k}).pCO2SW(MOORING.(omitmoors2{k}).month_since_1998==m),'omitnan');
        MOORING.(omitmoors2{k}).pCO2SW_mon_std(m)  = std(MOORING.(omitmoors2{k}).pCO2SW(MOORING.(omitmoors2{k}).month_since_1998==m),'omitnan');
    end
    % Plot climatology pCO2 and mooring measured pCO2
    figure; box on;
    set(gcf,'Position',[100, 400, 1680, 420])
    set(gca,'fontsize',20)
    title((omitmoor_names{k}),'fontsize',28);
    hold on
    %if strcmp('CB-06',omitmoors2{n}); lonidx = lonidx - 1; end
    %if strcmp('SEAK',omitmoors2{n}); latidx = latidx - 1; lonidx = lonidx - 2; end
    date = datenum([repmat(1998,276,1) SOCATv2021_grid.month_since_1998 ones(276,1)]);

    % borders of uncertainty envelopes
    t1=plot(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)) + ...
        U_tot_coast,'linewidth',1,'Color',clr1);
    b1=plot(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)) - ...
        U_tot_coast,'linewidth',1,'Color',clr1);
    t2=plot(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)) + ...
        U_tot_coast,'linewidth',1,'Color',clr4);
    b2=plot(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)) - ...
        U_tot_coast,'linewidth',1,'Color',clr4);
    t3=plot(date,MOORING.(omitmoors2{k}).pCO2SW_mon_mean + ...
        MOORING.(omitmoors2{k}).pCO2SW_mon_std,'linewidth',1,'Color','k');
    b3=plot(date,MOORING.(omitmoors2{k}).pCO2SW_mon_mean - ...
        MOORING.(omitmoors2{k}).pCO2SW_mon_std,'linewidth',1,'Color','k');
    t4=plot(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)) + 52.5,...
        'linewidth',1,'Color',clr2);
    b4=plot(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)) - 52.5,...
        'linewidth',1,'Color',clr2);
    % uncertainty envelopes
    f1=fill([date;flipud(date)],...
        [squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)) + U_tot_coast;...
         flipud(squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)) - U_tot_coast)],...
         clr1,'FaceAlpha',0.15,'LineStyle','none');
    f2=fill([date;flipud(date)],...
        [squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)) + U_tot_coast;...
         flipud(squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)) - U_tot_coast)],...
         clr4,'FaceAlpha',0.15,'LineStyle','none');
    % index to available mooring observations
    idx = find(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean) & ...
        ~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_std));
    a=diff(idx); % differences in index numbers
    b=find([a inf]>1); % ends of series
    c=diff([0 b]); % lengths of series
    d=cumsum(c); % ends of series
    idx1 = idx(1);
    for n=1:length(d)
        idx2 = idx(d(n));
        f3=fill([date(idx1:idx2);...
            flipud(date(idx1:idx2))],...
            [MOORING.(omitmoors2{k}).pCO2SW_mon_mean(idx1:idx2) + ...
             MOORING.(omitmoors2{k}).pCO2SW_mon_std((idx1:idx2)),...
             fliplr(MOORING.(omitmoors2{k}).pCO2SW_mon_mean(idx1:idx2) - ...
             MOORING.(omitmoors2{k}).pCO2SW_mon_std(idx1:idx2))],...
             'k','FaceAlpha',0.15,'LineStyle','none');
         if n < length(d)
             idx1 = idx(d(n)+1)
         end
    end
    f4=fill([LAR.date;flipud(LAR.date)],...
        [squeeze(LAR.pCO2(lonidx,latidx,:)) + 52.5;...
         flipud(squeeze(LAR.pCO2(lonidx,latidx,:)) - 52.5)],...
         clr2,'FaceAlpha',0.15,'LineStyle','none');
    % center lines
    p1=plot(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)),'linewidth',4,'Color',clr1);
    s1=scatter(date,squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:)),100,'o',...
        'MarkerFaceColor',clr1,'MarkerEdgeColor','none');
    p2=plot(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)),':','linewidth',4,'Color',clr4);
    s2=scatter(date,squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:)),100,'o',...
        'MarkerFaceColor',clr4,'MarkerEdgeColor','none');
    p3=plot(date,MOORING.(omitmoors2{k}).pCO2SW_mon_mean,'k','linewidth',4);
    s3=scatter(date,MOORING.(omitmoors2{k}).pCO2SW_mon_mean,100,'s',...
        'MarkerFaceColor','k','MarkerEdgeColor','none');
    p4=plot(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)),'Color',clr2,'linewidth',4);
    s4=scatter(LAR.date,squeeze(LAR.pCO2(lonidx,latidx,:)),100,'^',...
        'MarkerFaceColor',clr2,'MarkerEdgeColor','none');
    % dummy symbols for legends
    lg1 = plot(nan,nan, 'o-','linewidth',4,'MarkerSize',10,'MarkerEdgeColor','none',...
        'Color',clr1,'MarkerFaceColor',clr1,'DisplayName','RFR-CCS');
    lg2 = plot(nan,nan, 'o:','linewidth',4,'MarkerSize',10,'MarkerEdgeColor','none',...
        'Color',clr4,'MarkerFaceColor',clr4,'DisplayName','RFR-CCS-Eval');
    lg3 = plot(nan,nan, 's-','linewidth',4,'MarkerSize',10,'MarkerEdgeColor','none',...
        'Color','k','MarkerFaceColor','k','DisplayName','Mooring Observations');
    lg4 = plot(nan,nan, '^-','linewidth',4,'MarkerSize',10,'MarkerEdgeColor','none',...
        'Color',clr2,'MarkerFaceColor',clr2,'DisplayName','L17');
    % add legend
    legend([lg3 lg1 lg2 lg4],'location','northeast',...
                           'box','on','AutoUpdate','off');
    ylabel('{\itp}CO_{2(sw)} (\muatm)');
    xlabel('Year');
    datetick('x','yyyy');
    xlim([min(date(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))) ...
          max(date(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)))]);

    % Fit model between RF climatological pCO2 and mooring measurements
    RF=squeeze(SOCATv2021_grid.pco2_RF(lonidx,latidx,:));
    mdlRF = fitlm(MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)),...
        RF(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)));
    % Fit model between RF_val climatological pCO2 and mooring measurements
    RFval=squeeze(SOCATv2021_grid.pco2_RF_validate(lonidx,latidx,:));
    mdlRFval = fitlm(MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)),...
        RFval(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)));
    % Fit model between Laruelle climatological pCO2 and mooring measurements
    laruelle=squeeze(LAR.pCO2(lonidx,latidx,:));
    laruelle=[laruelle;nan(276-216,1)];
    mdlLAR = fitlm(MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean) & ~isnan(laruelle')),...
        laruelle(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean') & ~isnan(laruelle)));

    Mean_table(k,1) = SOCATv2021_grid.pco2_RF_annmean(lonidx,latidx);
    Mean_table(k,2) = SOCATv2021_grid.pco2_RF_validate_annmean(lonidx,latidx);
    Mean_table(k,3) = LAR.pCO2_annmean(lonidx,latidx);
    Mean_table(k,4) = mean(MOORING.(omitmoors2{k}).pCO2SW_monthly(:,1),'omitnan');

    Amplitude_table(k,1) = SOCATv2021_grid.pco2_RF_amp(lonidx,latidx);
    Amplitude_table(k,2) = SOCATv2021_grid.pco2_RF_validate_amp(lonidx,latidx);
    Amplitude_table(k,3) = LAR.pCO2_amp(lonidx,latidx);
    Amplitude_table(k,4) = MOORING.(omitmoors2{k}).pCO2SW_amplitude;

    R2_table(k,1) = mdlRF.Rsquared.Ordinary;
    R2_table(k,2) = mdlRFval.Rsquared.Ordinary;
    R2_table(k,3) = mdlLAR.Rsquared.Ordinary;

    Bias_table(k,1) = mean(RF(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))');
    Bias_table(k,2) = mean(RFval(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))');
    Bias_table(k,3) = mean(laruelle(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)' & ~isnan(laruelle))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)' & ~isnan(laruelle))');

    RMSE_table(k,1) = sqrt(sum((RF(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))').^2)./...
        sum(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)));
    RMSE_table(k,2) = sqrt(sum((RFval(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean))').^2)./...
        sum(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)));
    RMSE_table(k,3) = sqrt(sum((laruelle(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)' & ~isnan(laruelle))-...
        MOORING.(omitmoors2{k}).pCO2SW_mon_mean(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)' & ~isnan(laruelle))').^2)./...
        sum(~isnan(MOORING.(omitmoors2{k}).pCO2SW_mon_mean)' & ~isnan(laruelle)));

    exportgraphics(gcf,['Figures/mor_' omitmoors2{k} '.png']);
    
end

%% Plot map
% figure
% worldmap([latmin latmax],...
%     [lonmin lonmax]);
% %title('Amplitude of Seasonal pCO2 Cycle (Landschutzer merged climatology))');
% set(gcf,'Position',[617, 599, 700, 800])
% setm(gca,'ffacecolor',[0.94 0.97 1.0]);
% setm(gca,'fontsize',18)
% land = shaperead('landareas', 'UseGeoCoords', true);
% contourfm(SOCATv2021_grid.lat,SOCATv2021_grid.lon,SOCATv2021_grid.pco2_RF_amp,...
%     0:1:150,'LineStyle','none');
% geoshow(land, 'FaceColor', [0.2 0.2 0.2]);
% for n=1:numel(omitmoors2)
%     % Find closest point in climatology to mooring lat-lon
%     latidx = ...
%         find(abs(SOCATv2021_grid.lat(1,:)-nanmean(MOORING.(omitmoors2{n}).lat)) == ...
%         min(min(abs(SOCATv2021_grid.lat(1,:)-nanmean(MOORING.(omitmoors2{n}).lat)))));
%     lonidx = ...
%         find(abs(SOCATv2021_grid.lon(:,1)-nanmean(MOORING.(omitmoors2{n}).lon)) == ...
%         min(min(abs(SOCATv2021_grid.lon(:,1)-nanmean(MOORING.(omitmoors2{n}).lon)))));
%     %if strcmp('Southeast',moornames{n}); latidx = latidx - 1; end
%     scatterm(SOCATv2021_grid.lat(lonidx,latidx),SOCATv2021_grid.lon(lonidx,latidx),...
%         200,MOORING.(omitmoors2{n}).pCO2SW_amplitude,'filled',...
%         'MarkerEdgeColor','k','LineWidth',2);
%     if n==7
%     textm(SOCATv2021_grid.lat(lonidx,latidx)+1,SOCATv2021_grid.lon(lonidx,latidx)+3,...
%         omitmoors2{n},'Color','white','FontWeight','bold','fontsize',16);
%     else
%     textm(SOCATv2021_grid.lat(lonidx,latidx),SOCATv2021_grid.lon(lonidx,latidx)+3,...
%         omitmoors2{n},'Color','white','FontWeight','bold','fontsize',16);
%     end
% end
% c=colorbar; colormap(jet); caxis([0 150]);
% c.Label.String = 'pCO_{2} Seasonal Amplitude (\muatm)';
% c.FontSize = 18;