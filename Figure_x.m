
lat1 = 32:4:46;
lat2 = 34:4:48;

figure;
hold on;
set(gca,'fontsize',16);
%col = {'k' 'r' 'b' 'g'};

for n = 1:4

    %% Distance from shore index
    idx = SOCATv2020_grid.distance_from_shore(:,:,1) <= 100 & ...
          ~island(SOCATv2020_grid.lat,SOCATv2020_grid.lon) & ...
          SOCATv2020_grid.latitude(1,:,1) < lat2(n) & ...
          SOCATv2020_grid.latitude(1,:,1) > lat1(n);
    idx = repmat(idx,1,1,12);

    temp = SOCATv2020_grid.pco2_RF_clim;
    temp(~idx) = NaN;
    temp = squeeze(mean(mean(temp,1,'omitnan'),2,'omitnan'));
    matlab.lang.makeValidName(strcat('p',num2str(n))) = plot(1:12,temp,...
        'linewidth',2);
    
end

xlim([1 12]);
xlabel('Month of Year');
ylabel('pCO_{2} (\muatm)');
legend({['32-36 ' char(176) 'N'] ['36-40 ' char(176) 'N'] ...
        ['40-44 ' char(176) 'N'] ['44-48 ' char(176) 'N']},...
        'location','northwest');
hold off;

figure;
hold on;
set(gca,'fontsize',16);
%col = {'k' 'r' 'b' 'g'};

for n = 1:4

    %% Distance from shore index
    idx = SOCATv2020_grid.distance_from_shore(:,:,1) <= 100 & ...
          ~island(SOCATv2020_grid.lat,SOCATv2020_grid.lon) & ...
          SOCATv2020_grid.latitude(1,:,1) < lat2(n) & ...
          SOCATv2020_grid.latitude(1,:,1) > lat1(n);
    idx = repmat(idx,1,1,12);

    temp = SOCATv2020_grid.Fco2_RF_ERA5_avg;
    temp(~idx) = NaN;
    temp = squeeze(mean(mean(temp,1,'omitnan'),2,'omitnan'));
    matlab.lang.makeValidName(strcat('p',num2str(n))) = plot(1:12,temp,...
        'linewidth',2);
    
end

xlim([1 12]);
xlabel('Month of Year');
ylabel('FCO_{2} (mol m^{-2} yr^{-1})');
legend({['32-36 ' char(176) 'N'] ['36-40 ' char(176) 'N'] ...
        ['40-44 ' char(176) 'N'] ['44-48 ' char(176) 'N']},...
        'location','northwest');
hold off;
